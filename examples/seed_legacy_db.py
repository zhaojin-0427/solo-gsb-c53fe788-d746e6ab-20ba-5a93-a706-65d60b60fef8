#!/usr/bin/env python3
"""Seed an *old-version* database file for the upgrade walkthrough.

The historical service stored only ``requests``/``versions`` snapshots (no
``read_registry``).  Point this script at the data volume, start the current
image against the same file, and the service performs the crash-safe one-time
upgrade; the read ids below (``legacy-r0..legacy-r3``) then become
unrecoverable tombstones and reusing them returns HTTP 409.

Usage::

    python3 examples/seed_legacy_db.py ./data/umi.db
    # (inside the compose volume: UMI_DB_PATH=/data/umi.db)
"""
from __future__ import annotations

import json
import os
import sqlite3
import sys

SNAPSHOT = {
    "version": 1,
    "request_id": "legacy-run-001",
    "sample_id": "legacy-sample",
    "consensus": [
        {
            "locus": "IGH-A",
            "representative_umi": "AAAAA",
            "sequence": "ACGT",
            "depth": 4,
            "ambiguous": False,
        }
    ],
    "clusters": [
        {
            "locus": "IGH-A",
            "representative_umi": "AAAAA",
            "members": ["AAAAA"],
            "ambiguous": False,
            "depth": 4,
            "read_ids": [
                "legacy-r0",
                "legacy-r1",
                "legacy-r2",
                "legacy-r3",
            ],
        }
    ],
    "eliminated": [],
    "stats": {
        "loci": 1,
        "raw_reads": 4,
        "distinct_reads": 4,
        "duplicate_read_rows": 0,
        "clusters": 1,
        "consensus": 1,
        "eliminated": 0,
        "ambiguous_clusters": 0,
    },
}

OLD_SCHEMA = """
CREATE TABLE requests (
    request_id   TEXT PRIMARY KEY,
    sample_id    TEXT NOT NULL,
    payload_hash TEXT NOT NULL,
    version      INTEGER NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
);
CREATE TABLE versions (
    version      INTEGER PRIMARY KEY AUTOINCREMENT,
    request_id   TEXT NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
);
"""


def seed(path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    blob = json.dumps(SNAPSHOT, separators=(",", ":"), ensure_ascii=False)
    conn = sqlite3.connect(path)
    try:
        existing = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()
        if existing:
            raise SystemExit(
                f"refusing to seed {path!r}: file already contains tables "
                f"{sorted(r[0] for r in existing)}"
            )
        conn.executescript(OLD_SCHEMA)
        conn.execute(
            "INSERT INTO versions (version, request_id, created_at, result_json) "
            "VALUES (1, 'legacy-run-001', '2026-09-01T00:00:00Z', ?)",
            (blob,),
        )
        conn.execute(
            "INSERT INTO requests (request_id, sample_id, payload_hash, version, "
            "created_at, result_json) VALUES ('legacy-run-001', 'legacy-sample', "
            "'legacy-unknown-payload-hash', 1, '2026-09-01T00:00:00Z', ?)",
            (blob,),
        )
        conn.commit()
    finally:
        conn.close()
    print(f"seeded legacy database at {path} (version 1: legacy-run-001)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    seed(sys.argv[1])
