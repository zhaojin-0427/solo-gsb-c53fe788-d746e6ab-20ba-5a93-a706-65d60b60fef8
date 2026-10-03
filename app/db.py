"""SQLite persistence for immutable consensus versions.

Stored versions are append-only and never updated or deleted.  A
``BEGIN IMMEDIATE`` transaction plus a busy timeout makes same-``request_id``
idempotency *and* cross-request ``read_id`` content enforcement safe under
concurrent writers.

Tables
------
requests / versions
    One immutable snapshot per accepted request / per version number.

read_registry
    Persistent record of the canonical content (locus, umi, sequence and
    qualities hash) of every ``read_id`` ever accepted.  A ``read_id`` is
    global: submitting the same id from another ``request_id`` with different
    content is rejected with 409 inside the submit transaction.  Rows with
    ``content_hash IS NULL`` are *tombstones* migrated from a pre-upgrade
    database whose stored snapshots do not retain the original read payload;
    such ids can never be reused (see ``_migrate_legacy_reads``).
"""
from __future__ import annotations

import contextlib
import json
import os
import sqlite3
from typing import Dict, List, Optional, Tuple

from .consensus import build_consensus, payload_hash, read_content_hash
from .schemas import ConsensusRequest, ReadIn

DEFAULT_DB_PATH = os.environ.get("UMI_DB_PATH", "/data/umi.db")

# Marker used for created_at on rows back-filled from legacy snapshots: the
# values are not real submissions and must never be treated as content.
_LEGACY_MARK = "legacy-unknown"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS requests (
    request_id   TEXT PRIMARY KEY,
    sample_id    TEXT NOT NULL,
    payload_hash TEXT NOT NULL,
    version      INTEGER NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS versions (
    version      INTEGER PRIMARY KEY AUTOINCREMENT,
    request_id   TEXT NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS read_registry (
    read_id      TEXT PRIMARY KEY,
    locus        TEXT,
    umi          TEXT,
    sequence     TEXT,
    qualities    TEXT,
    content_hash TEXT,
    version      INTEGER NOT NULL,
    created_at   TEXT NOT NULL
);
"""


class ConflictError(Exception):
    """Base class for conflicts surfaced to the API as HTTP 409."""


class RequestPayloadConflict(ConflictError):
    """Same request_id replayed with a different payload (HTTP 409)."""


class ReadContentConflict(ConflictError):
    """A read_id was already persisted with different content (HTTP 409)."""


def _read_content_key(read: ReadIn) -> str:
    return read_content_hash(read)


class Database:
    def __init__(self, path: str = DEFAULT_DB_PATH) -> None:
        self.path = path
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)
        self._init_schema()

    @contextlib.contextmanager
    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path, timeout=30.0)
        try:
            conn.row_factory = sqlite3.Row
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA foreign_keys=ON")
            conn.execute("PRAGMA busy_timeout=30000")
            yield conn
        finally:
            conn.close()

    def _init_schema(self) -> None:
        with self._connect() as conn:
            # Detect a pre-upgrade database before creating the new table:
            # once read_registry exists the upgrade has already happened.
            legacy = conn.execute(
                "SELECT name FROM sqlite_master "
                "WHERE type='table' AND name='read_registry'"
            ).fetchone() is None
            conn.executescript(_SCHEMA)
            if legacy:
                self._migrate_legacy_reads(conn)
            conn.commit()

    @staticmethod
    def _migrate_legacy_reads(conn: sqlite3.Connection) -> None:
        """Back-fill read_ids found in pre-upgrade snapshots as tombstones.

        Old ``versions``/``requests`` rows only stored the computed result,
        whose ``clusters[].read_ids`` retain the ids but not the original
        locus/umi/sequence/qualities payload.  We must keep those versions
        readable unchanged, yet we also must not fabricate a content record
        for ids whose real content is unrecoverable: each id is therefore
        inserted with NULL content columns, which rejects every future reuse.
        """
        rows = conn.execute(
            "SELECT result_json FROM versions UNION ALL "
            "SELECT result_json FROM requests"
        ).fetchall()
        seen = set()
        for row in rows:
            try:
                snapshot = json.loads(row["result_json"])
            except (ValueError, TypeError):
                continue
            for cluster in snapshot.get("clusters", []):
                for read_id in cluster.get("read_ids", []):
                    if read_id not in seen:
                        seen.add(read_id)
                        conn.execute(
                            "INSERT OR IGNORE INTO read_registry "
                            "(read_id, locus, umi, sequence, qualities, "
                            "content_hash, version, created_at) "
                            "VALUES (?, NULL, NULL, NULL, NULL, NULL, ?, ?)",
                            (read_id, snapshot.get("version", 0), _LEGACY_MARK),
                        )

    # ------------------------------------------------------------------
    def get_by_request(self, request_id: str) -> Optional[dict]:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT result_json FROM requests WHERE request_id = ?",
                (request_id,),
            ).fetchone()
        return json.loads(row["result_json"]) if row else None

    def get_version(self, version: int) -> Optional[dict]:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT result_json FROM versions WHERE version = ?",
                (version,),
            ).fetchone()
        return json.loads(row["result_json"]) if row else None

    def latest_version(self) -> Optional[int]:
        with self._connect() as conn:
            row = conn.execute("SELECT MAX(version) AS v FROM versions").fetchone()
        return row["v"] if row and row["v"] is not None else None

    # ------------------------------------------------------------------
    @staticmethod
    def _unique_reads(req: ConsensusRequest) -> List[ReadIn]:
        """Collapse exact duplicate rows (same id *and* content)."""
        unique: Dict[Tuple, ReadIn] = {}
        for read in req.reads:
            unique.setdefault(
                (
                    read.read_id,
                    read.locus,
                    read.umi,
                    read.sequence,
                    tuple(read.qualities),
                ),
                read,
            )
        return list(unique.values())

    def submit(self, req: ConsensusRequest, created_at: str) -> Tuple[dict, bool]:
        """Insert a new version. Returns (result, replayed).

        All checks run atomically inside one ``BEGIN IMMEDIATE`` transaction:

        * same ``request_id`` + same payload hash -> replay stored version;
        * same ``request_id`` + different payload  -> RequestPayloadConflict;
        * a ``read_id`` already recorded (by any request) with different
          content -> ReadContentConflict; no version is allocated and no
          existing snapshot or registry row is modified;
        * a ``read_id`` matching a legacy tombstone (content unrecoverable
          from a pre-upgrade database) -> ReadContentConflict, never reused;
        * identical content under a new ``request_id`` is allowed and creates
          a new version.
        """
        phash = payload_hash(req)
        with self._connect() as conn:
            conn.execute("BEGIN IMMEDIATE")
            try:
                row = conn.execute(
                    "SELECT payload_hash, result_json FROM requests "
                    "WHERE request_id = ?",
                    (req.request_id,),
                ).fetchone()
                if row is not None:
                    conn.commit()
                    if row["payload_hash"] != phash:
                        raise RequestPayloadConflict(
                            f"request_id={req.request_id!r} already submitted with a "
                            "different payload"
                        )
                    return json.loads(row["result_json"]), True

                reads = self._unique_reads(req)

                # Cross-request read_id content check, against the durable
                # registry.  Exact duplicate rows were collapsed above; the
                # request-level validator (find_readid_conflicts) already
                # rejects same-id/different-content rows within one request.
                content_by_id = {r.read_id: r for r in reads}
                placeholders = ",".join("?" for _ in content_by_id)
                stored = conn.execute(
                    f"SELECT read_id, content_hash FROM read_registry "
                    f"WHERE read_id IN ({placeholders})",
                    tuple(content_by_id),
                ).fetchall()
                for existing in stored:
                    incoming = content_by_id[existing["read_id"]]
                    if existing["content_hash"] != _read_content_key(incoming):
                        raise ReadContentConflict(
                            f"read_id={existing['read_id']!r} was already submitted "
                            "with different content (locus/umi/sequence/qualities)"
                            if existing["content_hash"] is not None
                            else f"read_id={existing['read_id']!r} appeared in a "
                            "pre-upgrade result whose original content cannot be "
                            "recovered; the read_id cannot be reused"
                        )

                result = build_consensus(req)
                cur = conn.execute(
                    "INSERT INTO versions (request_id, created_at, result_json) "
                    "VALUES (?, ?, '{}')",
                    (req.request_id, created_at),
                )
                version = cur.lastrowid
                result["version"] = version
                result_blob = json.dumps(
                    result, separators=(",", ":"), ensure_ascii=False
                )
                conn.execute(
                    "UPDATE versions SET result_json = ? WHERE version = ?",
                    (result_blob, version),
                )
                conn.execute(
                    "INSERT INTO requests "
                    "(request_id, sample_id, payload_hash, version, created_at, "
                    "result_json) VALUES (?, ?, ?, ?, ?, ?)",
                    (
                        req.request_id,
                        req.sample_id,
                        phash,
                        version,
                        created_at,
                        result_blob,
                    ),
                )

                # Record each new read_id's content.  ids whose content
                # already matched are skipped (INSERT OR IGNORE keeps the
                # original first-seen row intact).
                for read in reads:
                    conn.execute(
                        "INSERT OR IGNORE INTO read_registry "
                        "(read_id, locus, umi, sequence, qualities, content_hash, "
                        "version, created_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        (
                            read.read_id,
                            read.locus,
                            read.umi,
                            read.sequence,
                            json.dumps(list(read.qualities), separators=(",", ":")),
                            _read_content_key(read),
                            version,
                            created_at,
                        ),
                    )
                conn.commit()
                return result, False
            except Exception:
                conn.rollback()
                raise
