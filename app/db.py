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

schema_migrations
    Completion markers for one-time upgrades.  The mere existence of a table
    created by an upgrade proves nothing: DDL is visible once committed while
    the data back-fill following it may still be in flight.  A migration is
    considered complete only when its marker row is committed **in the same
    write transaction** as its back-fill.  See ``_init_schema`` for the
    crash-safe / concurrent-startup protocol.
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

# Name of the one-time upgrade that adds read_registry and back-fills legacy
# read ids.  Applied == a committed row with this name in schema_migrations.
_READ_REGISTRY_MIGRATION = "0001_read_registry"

_SCHEMA_STATEMENTS = (
    """
CREATE TABLE IF NOT EXISTS requests (
    request_id   TEXT PRIMARY KEY,
    sample_id    TEXT NOT NULL,
    payload_hash TEXT NOT NULL,
    version      INTEGER NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
)
""",
    """
CREATE TABLE IF NOT EXISTS versions (
    version      INTEGER PRIMARY KEY AUTOINCREMENT,
    request_id   TEXT NOT NULL UNIQUE,
    created_at   TEXT NOT NULL,
    result_json  TEXT NOT NULL
)
""",
    """
CREATE TABLE IF NOT EXISTS read_registry (
    read_id      TEXT PRIMARY KEY,
    locus        TEXT,
    umi          TEXT,
    sequence     TEXT,
    qualities    TEXT,
    content_hash TEXT,
    version      INTEGER NOT NULL,
    created_at   TEXT NOT NULL
)
""",
    """
CREATE TABLE IF NOT EXISTS schema_migrations (
    name       TEXT PRIMARY KEY,
    applied_at TEXT NOT NULL
)
""",
)


def _utc_now_iso() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _is_busy_or_locked(exc: BaseException) -> bool:
    if not isinstance(exc, sqlite3.OperationalError):
        return False
    msg = str(exc).lower()
    return "database is locked" in msg or "database table is locked" in msg


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
        # Set once the upgrade marker is visible; submit() re-checks it before
        # opening its own write transaction as defence against any process
        # initialised against an incompletely upgraded file.
        self._migrations_complete = False
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

    def _migration_applied(self, conn: sqlite3.Connection, name: str) -> bool:
        exists = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' "
            "AND name='schema_migrations'"
        ).fetchone()
        if exists is None:
            return False
        return conn.execute(
            "SELECT 1 FROM schema_migrations WHERE name = ?", (name,)
        ).fetchone() is not None

    def _init_schema(self) -> None:
        """Create tables and (once) upgrade a pre-upgrade database.

        Crash-safe / concurrent-startup protocol:

        * Table creation alone is *not* treated as completion.  An upgrade is
          complete only when its back-fill *and* its ``schema_migrations``
          marker commit together in one write transaction; on failure the
          whole transaction (including ``CREATE TABLE``) is rolled back.
        * DDL is issued as plain ``execute`` calls, never ``executescript``
          (the latter force-commits any pending transaction first, which is
          exactly how a created-but-empty ``read_registry`` could survive a
          crash in the old code).
        * ``BEGIN IMMEDIATE`` takes the SQLite write lock, so simultaneous
          starts are serialised: the loser of the lock waits, then re-checks
          the marker and runs nothing if the winner finished the upgrade.
        """
        for _attempt in range(2):
            with self._connect() as conn:
                if self._migration_applied(conn, _READ_REGISTRY_MIGRATION):
                    # Fresh database or a fully upgraded one: just make sure
                    # every table exists.  CREATE TABLE IF NOT EXISTS inside a
                    # transaction commits atomically and changes nothing that
                    # is already there.
                    conn.execute("BEGIN IMMEDIATE")
                    try:
                        for statement in _SCHEMA_STATEMENTS:
                            conn.execute(statement)
                        conn.commit()
                    except Exception:
                        conn.rollback()
                        raise
                    self._migrations_complete = True
                    return

                # Upgrade needed (fresh DB or a pre-upgrade database).  The
                # IMMEDIATE lock makes concurrent initialisations wait here
                # until the first one commits or rolls back.
                conn.execute("BEGIN IMMEDIATE")
                try:
                    # Re-check under the lock: another process/thread may have
                    # completed the whole upgrade while we were waiting.
                    if self._migration_applied(
                        conn, _READ_REGISTRY_MIGRATION
                    ):
                        conn.commit()
                        self._migrations_complete = True
                        return
                    for statement in _SCHEMA_STATEMENTS:
                        conn.execute(statement)
                    self._migrate_legacy_reads(conn)
                    # Completion marker in the SAME transaction as the
                    # back-fill: a crash/failure before this commit leaves no
                    # marker, so the next start performs the upgrade again.
                    conn.execute(
                        "INSERT INTO schema_migrations (name, applied_at) "
                        "VALUES (?, ?)",
                        (_READ_REGISTRY_MIGRATION, _utc_now_iso()),
                    )
                    conn.commit()
                    self._migrations_complete = True
                    return
                except Exception as exc:
                    conn.rollback()
                    if not _is_busy_or_locked(exc):
                        raise
                    # SQLITE_BUSY despite busy_timeout: release this
                    # connection and retry; the holder will finish and the
                    # post-lock marker check makes the retry a no-op.
        raise RuntimeError("could not complete database migration (database locked)")

    @staticmethod
    def _migrate_legacy_reads(conn: sqlite3.Connection) -> None:
        """Back-fill read_ids found in pre-upgrade snapshots as tombstones.

        Old ``versions``/``requests`` rows only stored the computed result,
        whose ``clusters[].read_ids`` retain the ids but not the original
        locus/umi/sequence/qualities payload.  We must keep those versions
        readable unchanged, yet we also must not fabricate a content record
        for ids whose real content is unrecoverable: each id is therefore
        inserted with NULL content columns, which rejects every future reuse.

        INSERT OR IGNORE also guarantees a tombstone never overwrites a row
        carrying real content (e.g. a partially upgraded file created by the
        old, non-atomic code that already holds new submissions).
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
        # A submit must never slip past an unfinished upgrade.  Taking the
        # write lock here queues behind a migration still running in another
        # process (busy_timeout); once the lock is acquired the marker is
        # re-checked, so a completed upgrade lets the submit proceed while an
        # incomplete one aborts it instead of registering reads against a
        # half-migrated registry.
        if not self._migrations_complete:
            with self._connect() as guard:
                guard.execute("BEGIN IMMEDIATE")
                try:
                    applied = self._migration_applied(
                        guard, _READ_REGISTRY_MIGRATION
                    )
                    guard.commit()
                except Exception:
                    guard.rollback()
                    raise
                if not applied:
                    raise RuntimeError("database migration has not completed")
            self._migrations_complete = True

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
