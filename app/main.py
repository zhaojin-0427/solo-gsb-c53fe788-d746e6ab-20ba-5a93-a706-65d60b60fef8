"""FastAPI application exposing the UMI consensus denoising API."""
from __future__ import annotations

import os
from contextlib import asynccontextmanager
from datetime import datetime, timezone

from fastapi import Depends, FastAPI, HTTPException, Request, Response, status
from fastapi.responses import JSONResponse

from .consensus import find_readid_conflicts
from .db import (
    DEFAULT_DB_PATH,
    Database,
    ReadContentConflict,
    RequestPayloadConflict,
)
from .schemas import ConsensusRequest, ConsensusResult


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Read the path at startup so tests/containers can point UMI_DB_PATH
    # at an alternative location.
    app.state.db = Database(os.environ.get("UMI_DB_PATH", DEFAULT_DB_PATH))
    yield


app = FastAPI(
    title="UMI Molecular Consensus Denoising API",
    version="1.1.1",
    lifespan=lifespan,
    description=(
        "Cluster reads by UMI (frequency/lexicographic ordering, Hamming "
        "distance <= 1 assignment, multi-hit clusters ambiguous) and build a "
        "Phred-weighted per-position consensus. Versions are immutable. "
        "read_id content is recorded durably: the same id may never carry "
        "different content, even across different request_ids."
    ),
)


def get_db(request: Request) -> Database:
    return request.app.state.db


@app.exception_handler(ValueError)
async def _value_error_handler(_: Request, exc: ValueError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"error": str(exc)})


@app.get("/healthz")
def healthz(db: Database = Depends(get_db)) -> dict:
    latest = db.latest_version()
    return {"status": "ok", "latest_version": latest}


@app.post(
    "/consensus",
    response_model=ConsensusResult,
    status_code=status.HTTP_200_OK,
    tags=["consensus"],
)
def post_consensus(
    req: ConsensusRequest,
    response: Response,
    db: Database = Depends(get_db),
) -> dict:
    """Create (or idempotently replay) a versioned consensus computation."""
    conflicts = find_readid_conflicts(req)
    if conflicts:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail={
                "error": "read_id content conflict",
                "conflicts": conflicts,
            },
        )

    try:
        result, replayed = db.submit(req, utc_now_iso())
    except ReadContentConflict as exc:
        # Cross-request: a read_id already persisted with different content,
        # or a read_id carried over from a pre-upgrade database whose
        # original payload cannot be recovered.
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail={"error": str(exc)},
        )
    except RequestPayloadConflict as exc:
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail={"error": str(exc), "request_id": req.request_id},
        )
    response.headers["X-Idempotent-Replay"] = "true" if replayed else "false"
    response.headers["X-Version"] = str(result["version"])
    return result


@app.get(
    "/consensus/{version}",
    response_model=ConsensusResult,
    tags=["consensus"],
)
def get_consensus(version: int, db: Database = Depends(get_db)) -> dict:
    """Return the immutable snapshot stored under ``version``."""
    result = db.get_version(version)
    if result is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail={"error": f"version {version} not found"},
        )
    return result
