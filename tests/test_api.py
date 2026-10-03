"""End-to-end API tests: idempotency, conflicts, immutability."""
import json
import sqlite3
import threading
import time

import pytest
from fastapi.testclient import TestClient

PAYLOAD = {
    "request_id": "req-1",
    "sample_id": "sample-A",
    "reads": [
        {
            "locus": "IGH-A",
            "read_id": f"r{i}",
            "umi": "AAAAA",
            "sequence": "ACGT",
            "qualities": [30, 30, 30, 30],
        }
        for i in range(4)
    ],
}


def test_post_creates_version_and_get_is_identical(client: TestClient):
    resp = client.post("/consensus", json=PAYLOAD)
    assert resp.status_code == 200, resp.text
    assert resp.headers["X-Idempotent-Replay"] == "false"
    body = resp.json()
    assert body["version"] == 1
    assert body["consensus"][0]["sequence"] == "ACGT"
    assert body["request_id"] == "req-1"

    got = client.get(f"/consensus/{body['version']}")
    assert got.status_code == 200
    assert got.json() == body


def test_same_request_same_payload_is_idempotent(client: TestClient):
    first = client.post("/consensus", json=PAYLOAD)
    assert first.status_code == 200

    # Reorder reads: still the same canonical payload.
    reordered = dict(PAYLOAD, reads=list(reversed(PAYLOAD["reads"])))
    second = client.post("/consensus", json=reordered)
    assert second.status_code == 200
    assert second.headers["X-Idempotent-Replay"] == "true"
    assert second.json() == first.json()

    assert client.get("/consensus/2").status_code == 404  # no version 2
    assert client.get("/healthz").json()["latest_version"] == 1


def test_same_request_different_payload_is_409(client: TestClient):
    assert client.post("/consensus", json=PAYLOAD).status_code == 200
    changed = {
        **PAYLOAD,
        "reads": [
            {**PAYLOAD["reads"][0], "sequence": "ACGG"},
            *PAYLOAD["reads"][1:],
        ],
    }
    resp = client.post("/consensus", json=changed)
    assert resp.status_code == 409
    assert "different payload" in resp.json()["detail"]["error"]


def test_readid_content_conflict_is_409(client: TestClient):
    payload = {
        "request_id": "req-2",
        "sample_id": "s",
        "reads": [
            {
                "locus": "L1",
                "read_id": "dup",
                "umi": "AAAAA",
                "sequence": "ACGT",
                "qualities": [30] * 4,
            },
            {
                "locus": "L1",
                "read_id": "dup",
                "umi": "AAAAT",
                "sequence": "ACGT",
                "qualities": [30] * 4,
            },
            {
                "locus": "L1",
                "read_id": "other",
                "umi": "AAAAA",
                "sequence": "ACGT",
                "qualities": [30] * 4,
            },
        ],
    }
    resp = client.post("/consensus", json=payload)
    assert resp.status_code == 409
    assert resp.json()["detail"]["conflicts"][0]["read_id"] == "dup"


def test_identical_duplicate_rows_are_allowed(client: TestClient):
    read = {
        "locus": "L1",
        "read_id": "r1",
        "umi": "AAAAA",
        "sequence": "ACGT",
        "qualities": [30] * 4,
    }
    payload = {
        "request_id": "req-3",
        "sample_id": "s",
        "reads": [dict(read), dict(read)],
    }
    resp = client.post("/consensus", json=payload)
    assert resp.status_code == 200
    assert resp.json()["stats"]["duplicate_read_rows"] == 1


def test_validation_errors_are_422(client: TestClient):
    bad = {
        "request_id": "req-4",
        "sample_id": "s",
        "reads": [
            {
                "locus": "L1",
                "read_id": "r1",
                "umi": "AAAAA",
                "sequence": "ACGT",
                "qualities": [30, 30],  # length mismatch
            }
        ],
    }
    assert client.post("/consensus", json=bad).status_code == 422

    bad_q = {
        "request_id": "req-5",
        "sample_id": "s",
        "reads": [
            {
                "locus": "L1",
                "read_id": "r1",
                "umi": "AAAAA",
                "sequence": "ACGT",
                "qualities": [30, 50, 30, 30],  # out of range
            }
        ],
    }
    assert client.post("/consensus", json=bad_q).status_code == 422


def test_versions_are_immutable_and_independent(client: TestClient):
    v1 = client.post("/consensus", json=PAYLOAD).json()
    second_payload = {
        "request_id": "req-2",
        "sample_id": "sample-B",
        "reads": [
            {
                "locus": "IGH-B",
                "read_id": f"b{i}",
                "umi": "TTTTT",
                "sequence": "TTTT",
                "qualities": [40] * 4,
            }
            for i in range(3)
        ],
    }
    v2 = client.post("/consensus", json=second_payload).json()
    assert v2["version"] == 2

    # Re-reading v1 must still return the original snapshot.
    assert client.get("/consensus/1").json() == v1
    assert client.get("/consensus/2").json() == v2
    assert client.get("/consensus/999").status_code == 404


# --------------------------------------------------------------------------
# Cross-request read_id content enforcement
# --------------------------------------------------------------------------

def _read(read_id, *, locus="L1", umi="AAAAA", sequence="ACGT", quals=None):
    return {
        "locus": locus,
        "read_id": read_id,
        "umi": umi,
        "sequence": sequence,
        "qualities": [30] * len(sequence) if quals is None else quals,
    }


def test_same_content_under_new_request_creates_new_version(client: TestClient):
    first = {
        "request_id": "cross-req-1",
        "sample_id": "s",
        "reads": [_read("shared-1"), _read("shared-2"), _read("shared-3"),
                  _read("shared-4")],
    }
    resp1 = client.post("/consensus", json=first)
    assert resp1.status_code == 200, resp1.text
    assert resp1.headers["X-Idempotent-Replay"] == "false"

    # A different request_id reusing the exact same read content is allowed
    # and creates its own version.
    second = dict(first, request_id="cross-req-2")
    resp2 = client.post("/consensus", json=second)
    assert resp2.status_code == 200, resp2.text
    assert resp2.headers["X-Idempotent-Replay"] == "false"
    assert resp2.json()["version"] == 2
    assert client.get("/healthz").json()["latest_version"] == 2


@pytest.mark.parametrize(
    "changed",
    [
        {"locus": "L2"},
        {"umi": "AAAAT"},
        {"sequence": "ACGG"},
        {"qualities": [30, 30, 30, 25]},
    ],
)
def test_cross_request_readid_conflict_is_409(client: TestClient, changed):
    first = {
        "request_id": "cross-A",
        "sample_id": "s",
        "reads": [_read(f"r{i}") for i in range(4)],
    }
    assert client.post("/consensus", json=first).status_code == 200

    second_reads = [_read(f"r{i}") for i in range(4)]
    second_reads[0] = {**second_reads[0], **changed}
    second = {"request_id": "cross-B", "sample_id": "s", "reads": second_reads}
    resp = client.post("/consensus", json=second)
    assert resp.status_code == 409
    assert "already submitted" in resp.json()["detail"]["error"]

    # No new version, no snapshot mutation.
    assert client.get("/healthz").json()["latest_version"] == 1
    snap = client.get("/consensus/1").json()
    assert sorted(snap["clusters"][0]["read_ids"]) == ["r0", "r1", "r2", "r3"]

    # The conflicting request_id itself was never registered, so retrying it
    # with the original content is a fresh successful submit.
    fixed = {"request_id": "cross-B", "sample_id": "s",
             "reads": [_read(f"r{i}") for i in range(4)]}
    retry = client.post("/consensus", json=fixed)
    assert retry.status_code == 200
    assert retry.json()["version"] == 2


def test_cross_request_conflict_does_not_poison_other_reads(
    client: TestClient, tmp_path
):
    first = {
        "request_id": "cross-A",
        "sample_id": "s",
        "reads": [_read(f"r{i}") for i in range(4)],
    }
    client.post("/consensus", json=first)

    # r0 conflicts; brand-new r9/r10 must NOT be registered by the failed tx.
    bad_reads = [
        _read("r0", sequence="ACGG"),
        _read("r9"),
        _read("r10"),
        _read("r11"),
    ]
    bad = {"request_id": "cross-B", "sample_id": "s", "reads": bad_reads}
    assert client.post("/consensus", json=bad).status_code == 409

    # A third request that reuses the "new" ids with *different* content must
    # still succeed in defining them: the rolled-back tx left no trace.
    third_reads = [
        _read("r9", locus="L2"),
        _read("r10", locus="L2"),
        _read("r11", locus="L2"),
        _read("r12", locus="L2"),
    ]
    third = {"request_id": "cross-C", "sample_id": "s", "reads": third_reads}
    resp = client.post("/consensus", json=third)
    assert resp.status_code == 200, resp.text
    assert resp.json()["version"] == 2


def test_concurrent_conflicting_submits_only_one_wins(client: TestClient):
    barrier = threading.Barrier(2)
    outcomes = {}

    def submit(thread_name, sequence):
        payload = {
            "request_id": f"conc-{thread_name}",
            "sample_id": "s",
            "reads": [
                {
                    "locus": "L1",
                    "read_id": f"shared-{i}",
                    "umi": "AAAAA",
                    "sequence": sequence,
                    "qualities": [30] * 4,
                }
                for i in range(4)
            ],
        }
        barrier.wait()
        resp = client.post("/consensus", json=payload)
        outcomes[thread_name] = (resp.status_code, resp.text)

    t1 = threading.Thread(target=submit, args=("one", "ACGT"))
    t2 = threading.Thread(target=submit, args=("two", "ACGG"))
    t1.start()
    t2.start()
    t1.join()
    t2.join()

    statuses = sorted(code for code, _ in outcomes.values())
    assert statuses == [200, 409], outcomes
    assert client.get("/healthz").json()["latest_version"] == 1

    # Only the winning content is recorded: the loser's id reuse with the
    # winning content now succeeds (new request, identical content).
    winning = next(text for code, text in outcomes.values() if code == 200)
    winner_seq = json.loads(winning)["consensus"][0]["sequence"]
    followup = {
        "request_id": "conc-followup",
        "sample_id": "s",
        "reads": [
            {
                "locus": "L1",
                "read_id": f"shared-{i}",
                "umi": "AAAAA",
                "sequence": winner_seq,
                "qualities": [30] * 4,
            }
            for i in range(4)
        ],
    }
    ok = client.post("/consensus", json=followup)
    assert ok.status_code == 200, ok.text
    assert ok.json()["version"] == 2


# --------------------------------------------------------------------------
# Upgrade of a pre-upgrade database (no read_registry table)
# --------------------------------------------------------------------------

def _seed_legacy_db(path) -> dict:
    """Create a database in the old schema with one stored snapshot."""
    snapshot = {
        "version": 1,
        "request_id": "legacy-req",
        "sample_id": "legacy-sample",
        "consensus": [
            {
                "locus": "L1",
                "representative_umi": "AAAAA",
                "sequence": "ACGT",
                "depth": 4,
                "ambiguous": False,
            }
        ],
        "clusters": [
            {
                "locus": "L1",
                "representative_umi": "AAAAA",
                "members": ["AAAAA"],
                "ambiguous": False,
                "depth": 4,
                "read_ids": ["old-0", "old-1", "old-2", "old-3"],
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
    blob = json.dumps(snapshot, separators=(",", ":"))
    conn = sqlite3.connect(path)
    try:
        conn.executescript(
            """
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
        )
        conn.execute(
            "INSERT INTO versions (version, request_id, created_at, result_json) "
            "VALUES (1, 'legacy-req', '2026-09-01T00:00:00Z', ?)",
            (blob,),
        )
        conn.execute(
            "INSERT INTO requests (request_id, sample_id, payload_hash, version, "
            "created_at, result_json) VALUES ('legacy-req', 'legacy-sample', "
            "'oldhash', 1, '2026-09-01T00:00:00Z', ?)",
            (blob,),
        )
        conn.commit()
    finally:
        conn.close()
    return snapshot


def test_legacy_db_keeps_old_versions_readable_and_tombstones_read_ids(
    tmp_path, monkeypatch
):
    db_path = tmp_path / "legacy.db"
    legacy_snapshot = _seed_legacy_db(db_path)

    monkeypatch.setenv("UMI_DB_PATH", str(db_path))
    from app.main import app

    with TestClient(app) as c:
        # Original version is still readable byte-for-byte.
        got = c.get("/consensus/1")
        assert got.status_code == 200
        assert got.json() == legacy_snapshot
        assert c.get("/healthz").json()["latest_version"] == 1

        # Reusing a legacy read_id with any (even plausible) content is 409:
        # its original content cannot be reconstructed from the snapshot.
        reuse = {
            "request_id": "after-upgrade-1",
            "sample_id": "s",
            "reads": [_read(f"old-{i}") for i in range(4)],
        }
        resp = c.post("/consensus", json=reuse)
        assert resp.status_code == 409
        assert "pre-upgrade result" in resp.json()["detail"]["error"]
        assert c.get("/healthz").json()["latest_version"] == 1

        # Brand-new read_ids work normally and allocate a fresh version.
        fresh = {
            "request_id": "after-upgrade-2",
            "sample_id": "s",
            "reads": [_read(f"new-{i}") for i in range(4)],
        }
        ok = c.post("/consensus", json=fresh)
        assert ok.status_code == 200
        assert ok.json()["version"] == 2

        # ...and from there on cross-request enforcement applies normally.
        conflict = {
            "request_id": "after-upgrade-3",
            "sample_id": "s",
            "reads": [
                {**_read("new-0"), "sequence": "ACGG"},
                _read("new-1"),
                _read("new-2"),
                _read("new-3"),
            ],
        }
        assert c.post("/consensus", json=conflict).status_code == 409

        # Old snapshot remains untouched.
        assert c.get("/consensus/1").json() == legacy_snapshot


# --------------------------------------------------------------------------
# Recovery from an upgrade interrupted after read_registry was created but
# before the legacy back-fill committed (the old code treated the table's
# existence as completion and skipped the back-fill on restart).
# --------------------------------------------------------------------------

def _simulate_interrupted_upgrade(path) -> dict:
    """Legacy DB plus a committed, EMPTY read_registry and no marker."""
    snapshot = _seed_legacy_db(path)
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            """
            CREATE TABLE read_registry (
                read_id      TEXT PRIMARY KEY,
                locus        TEXT,
                umi          TEXT,
                sequence     TEXT,
                qualities    TEXT,
                content_hash TEXT,
                version      INTEGER NOT NULL,
                created_at   TEXT NOT NULL
            )
            """
        )
        conn.commit()
    finally:
        conn.close()
    return snapshot


def test_interrupted_upgrade_is_resumed_on_restart(tmp_path, monkeypatch):
    db_path = tmp_path / "interrupted.db"
    legacy_snapshot = _simulate_interrupted_upgrade(db_path)

    # The empty table must not count as "migrated": tombstones are missing.
    conn = sqlite3.connect(db_path)
    try:
        assert conn.execute("SELECT count(*) FROM read_registry").fetchone()[0] == 0
    finally:
        conn.close()

    monkeypatch.setenv("UMI_DB_PATH", str(db_path))
    from app.main import app

    with TestClient(app) as c:
        # Back-fill resumed: old snapshot intact and old ids are tombstones.
        assert c.get("/consensus/1").json() == legacy_snapshot

        reuse = {
            "request_id": "after-recovery-1",
            "sample_id": "s",
            "reads": [_read(f"old-{i}") for i in range(4)],
        }
        resp = c.post("/consensus", json=reuse)
        assert resp.status_code == 409
        assert "pre-upgrade result" in resp.json()["detail"]["error"]
        assert c.get("/healthz").json()["latest_version"] == 1

        # Brand-new ids still submit normally.
        fresh = {
            "request_id": "after-recovery-2",
            "sample_id": "s",
            "reads": [_read(f"new-{i}") for i in range(4)],
        }
        ok = c.post("/consensus", json=fresh)
        assert ok.status_code == 200
        assert ok.headers["X-Idempotent-Replay"] == "false"
        assert ok.json()["version"] == 2

        # Same request_id + same payload replays v2; different payload 409.
        replay = c.post("/consensus", json=fresh)
        assert replay.status_code == 200
        assert replay.headers["X-Idempotent-Replay"] == "true"
        assert replay.json()["version"] == 2

        changed = {
            **fresh,
            "reads": [
                {**fresh["reads"][0], "sequence": "ACGG"},
                *fresh["reads"][1:],
            ],
        }
        assert c.post("/consensus", json=changed).status_code == 409
        assert c.get("/healthz").json()["latest_version"] == 2

        assert c.get("/consensus/1").json() == legacy_snapshot

    # Restart against the same file: migration is idempotent and the marker
    # prevents any re-run; reuse still 409.
    with TestClient(app) as c2:
        assert c2.get("/consensus/2").status_code == 200
        again = c2.post(
            "/consensus",
            json={
                "request_id": "after-recovery-3",
                "sample_id": "s",
                "reads": [_read("old-0"), _read("old-1"),
                          _read("old-2"), _read("old-3")],
            },
        )
        assert again.status_code == 409
        assert c2.get("/healthz").json()["latest_version"] == 2


def test_migration_backfill_and_marker_commit_together(
    tmp_path, monkeypatch
):
    """A failing back-fill commits nothing: table and marker roll back."""
    from app import db as db_module

    db_path = tmp_path / "rollback.db"
    _seed_legacy_db(db_path)

    original = db_module.Database._migrate_legacy_reads

    def fail_after_backfill(conn):
        original(conn)
        raise RuntimeError("simulated crash before commit")

    monkeypatch.setattr(
        db_module.Database,
        "_migrate_legacy_reads",
        staticmethod(fail_after_backfill),
    )
    with pytest.raises(RuntimeError, match="simulated crash"):
        db_module.Database(str(db_path))

    # The entire migration transaction (DDL included) rolled back.
    conn = sqlite3.connect(db_path)
    try:
        tables = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        assert "read_registry" not in tables
        assert "schema_migrations" not in tables
    finally:
        conn.close()

    # A later start retries and completes the upgrade normally.
    monkeypatch.undo()
    db = db_module.Database(str(db_path))
    assert db._migrations_complete is True
    conn = sqlite3.connect(db_path)
    try:
        assert conn.execute(
            "SELECT name FROM schema_migrations"
        ).fetchall() == [("0001_read_registry",)]
        assert conn.execute(
            "SELECT count(*) FROM read_registry WHERE content_hash IS NULL"
        ).fetchone()[0] == 4
    finally:
        conn.close()


def test_concurrent_starts_serialise_migration(tmp_path, monkeypatch):
    """Two processes/threads upgrading at once: one winner, no double work."""
    from app import db as db_module

    db_path = tmp_path / "concurrent.db"
    _seed_legacy_db(db_path)

    original = db_module.Database._migrate_legacy_reads
    hold = threading.Event()

    def slow_backfill(conn):
        hold.wait(10)
        time.sleep(0.5)  # keep the write lock so the waiter is definitely queued
        original(conn)

    monkeypatch.setattr(
        db_module.Database,
        "_migrate_legacy_reads",
        staticmethod(slow_backfill),
    )

    outcomes = {}

    def init(name):
        try:
            db_module.Database(str(db_path))
            outcomes[name] = "ok"
        except Exception as exc:  # pragma: no cover - failure path
            outcomes[name] = repr(exc)

    t1 = threading.Thread(target=init, args=("first",))
    t2 = threading.Thread(target=init, args=("second",))
    t1.start()
    # Winner takes the BEGIN IMMEDIATE lock and enters the slowed back-fill.
    time.sleep(0.5)
    t2.start()
    # The waiter is now queued behind the winner's lock; release the winner.
    time.sleep(0.5)
    hold.set()
    t1.join(30)
    t2.join(30)

    assert outcomes == {"first": "ok", "second": "ok"}, outcomes
    conn = sqlite3.connect(db_path)
    try:
        assert conn.execute(
            "SELECT count(*) FROM schema_migrations"
        ).fetchone()[0] == 1
        assert conn.execute(
            "SELECT count(*) FROM read_registry"
        ).fetchone()[0] == 4
    finally:
        conn.close()


def test_tombstone_backfill_never_overwrites_real_content(tmp_path):
    """A registry row carrying real content survives a repeated back-fill."""
    from app import db as db_module

    db_path = tmp_path / "mixed.db"
    _seed_legacy_db(db_path)
    db = db_module.Database(str(db_path))
    assert db._migrations_complete is True

    # Re-running the back-fill directly (as a restart/resume would against a
    # file with mixed state) must not replace real rows with NULL tombstones.
    with db._connect() as conn:
        conn.execute("BEGIN IMMEDIATE")
        conn.execute(
            "INSERT OR REPLACE INTO read_registry "
            "(read_id, locus, umi, sequence, qualities, content_hash, "
            "version, created_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            ("old-0", "L1", "AAAAA", "ACGT", "[30,30,30,30]",
             "realhash", 2, "2026-10-03T00:00:00Z"),
        )
        db_module.Database._migrate_legacy_reads(conn)
        conn.commit()

    with db._connect() as conn:
        row = conn.execute(
            "SELECT locus, content_hash FROM read_registry WHERE read_id='old-0'"
        ).fetchone()
    assert row["locus"] == "L1"
    assert row["content_hash"] == "realhash"
