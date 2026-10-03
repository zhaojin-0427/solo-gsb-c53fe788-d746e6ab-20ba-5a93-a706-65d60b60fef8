"""UMI clustering and Phred-weighted consensus voting.

Pipeline for one request payload:

1. De-duplicate fully identical read rows; reject same ``read_id`` with
   different content as a client conflict (409).
2. Group reads by locus.
3. Within each locus, sort distinct UMIs by frequency (descending) then
   lexicographically.  Each UMI starts a new cluster when no existing seed
   representative is within Hamming distance 1.  A UMI hitting seeds from
   multiple clusters merges them and marks the merged cluster *ambiguous*.
4. Within a cluster, each position is a Phred-weighted vote
   ``weight(base) = sum 10**(Q/10)``.  A base is called only when the
   winning weight is at least twice the sum of every other base's weight
   *and* the position is covered by >= 3 reads; otherwise the position is
   ``N`` and that cluster's consensus is eliminated.
"""
from __future__ import annotations

import json
from collections import Counter, OrderedDict, defaultdict
from typing import Any, Dict, List, Tuple

from .schemas import ConsensusResult, ReadIn

BASES = ("A", "C", "G", "T")
MIN_COVERAGE = 3
WEIGHT_MARGIN = 2.0


# ---------------------------------------------------------------------------
# Payload normalisation / idempotency
# ---------------------------------------------------------------------------

def _read_tuple(r: ReadIn) -> Tuple[str, str, str, str, Tuple[int, ...]]:
    return (r.locus, r.read_id, r.umi, r.sequence, tuple(r.qualities))


def canonical_payload(req: ConsensusRequest) -> dict:
    """Order-independent canonical representation of a payload.

    Read order in the JSON array must not affect identity, so rows are
    sorted and duplicates collapsed.
    """
    rows = sorted(_read_tuple(r) for r in req.reads)
    return {
        "request_id": req.request_id,
        "sample_id": req.sample_id,
        "reads": [
            {
                "locus": a,
                "read_id": b,
                "umi": c,
                "sequence": d,
                "qualities": list(e),
            }
            for a, b, c, d, e in rows
        ],
    }


def payload_hash(req: ConsensusRequest) -> str:
    import hashlib

    blob = json.dumps(canonical_payload(req), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def read_content_hash(r: ReadIn) -> str:
    """Stable fingerprint of one read's content (everything but read_id).

    Two rows sharing a ``read_id`` are considered the same read only when
    locus, umi, sequence and qualities are byte-for-byte identical.
    """
    import hashlib

    blob = json.dumps(
        [r.locus, r.umi, r.sequence, list(r.qualities)],
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def find_readid_conflicts(req: ConsensusRequest) -> List[dict]:
    """Return rows describing same read_id carrying different content."""
    seen: Dict[str, Tuple[str, str, str, Tuple[int, ...]]] = {}
    conflicts: List[dict] = []
    for r in req.reads:
        sig = (r.locus, r.umi, r.sequence, tuple(r.qualities))
        if r.read_id in seen and seen[r.read_id] != sig:
            conflicts.append(
                {"read_id": r.read_id, "reason": "same read_id with different content"}
            )
        else:
            seen[r.read_id] = sig
    return conflicts


def dedupe_reads(req: ConsensusRequest) -> Tuple[List[ReadIn], int]:
    """Drop exact duplicate rows; return unique reads and dropped count."""
    unique: "OrderedDict[Tuple, ReadIn]" = OrderedDict()
    for r in req.reads:
        unique.setdefault(_read_tuple(r), r)
    dropped = len(req.reads) - len(unique)
    return list(unique.values()), dropped


# ---------------------------------------------------------------------------
# UMI clustering
# ---------------------------------------------------------------------------

def hamming_le1(a: str, b: str) -> bool:
    """True iff Hamming distance(a, b) <= 1. Different lengths never match."""
    if len(a) != len(b):
        return False
    diff = 0
    for x, y in zip(a, b):
        if x != y:
            diff += 1
            if diff > 1:
                return False
    return True


class _Cluster:
    __slots__ = ("seeds", "members", "ambiguous", "representative")

    def __init__(self, seed: str) -> None:
        self.seeds: List[str] = [seed]
        self.members: List[str] = [seed]
        self.ambiguous: bool = False
        self.representative: str = seed

    def merge(self, other: "_Cluster") -> None:
        if other is self:
            return
        # representative = highest-priority seed (the one created first)
        self.seeds.extend(other.seeds)
        for m in other.members:
            if m not in self.members:
                self.members.append(m)
        self.ambiguous = True


def cluster_locus(reads: List[ReadIn]) -> List[_Cluster]:
    """Cluster reads by UMI within one locus.

    Two phases, both driven by the frequency-desc / lexicographic-asc order:

    1. Select representatives ("seeds"): a UMI becomes a representative iff
       it is more than Hamming distance 1 away from every previously chosen
       representative.
    2. Assign every non-representative UMI to *all* representatives within
       Hamming distance 1.  Hitting a single representative joins that
       cluster; hitting two or more merges them and marks the whole merged
       cluster ambiguous.
    """
    umi_freq = Counter(r.umi for r in reads)
    ordered = sorted(umi_freq, key=lambda u: (-umi_freq[u], u))

    # Phase 1: representatives.
    seeds: List[str] = []
    non_seeds: List[str] = []
    for umi in ordered:
        if any(hamming_le1(umi, s) for s in seeds):
            non_seeds.append(umi)
        else:
            seeds.append(umi)

    clusters = {seed: _Cluster(seed) for seed in seeds}

    # Phase 2: assign non-seeds against the full representative set.
    for umi in non_seeds:
        # Deduplicate by identity: already-merged seeds alias the same object.
        hits: List[_Cluster] = []
        for s in seeds:
            if hamming_le1(umi, s):
                cl = clusters[s]
                if cl not in hits:
                    hits.append(cl)
        if len(hits) == 1:
            hits[0].members.append(umi)
            continue
        # len(hits) >= 2: the UMI hits multiple representatives -> merge
        # every hit cluster and flag the merged whole cluster ambiguous.
        target = hits[0]
        target.members.append(umi)
        for cl in hits[1:]:
            target.merge(cl)
        for seed in seeds:
            if clusters[seed] in hits[1:]:
                clusters[seed] = target  # alias to the merged cluster

    # Unique clusters in seed priority order (highest-priority seed first).
    result: List[_Cluster] = []
    for seed in seeds:
        cl = clusters[seed]
        if cl not in result:
            result.append(cl)
    return result


# ---------------------------------------------------------------------------
# Weighted consensus
# ---------------------------------------------------------------------------

def _vote_position(
    reads: List[ReadIn], pos: int
) -> Tuple[str, float, int, Dict[str, float], bool]:
    """Return (base or 'N', top_weight, coverage, weights_by_base, tie)."""
    weights: Dict[str, float] = {b: 0.0 for b in BASES}
    coverage = 0
    for r in reads:
        if pos < len(r.sequence):
            coverage += 1
            weights[r.sequence[pos]] += 10.0 ** (r.qualities[pos] / 10.0)

    ranked = sorted(weights.items(), key=lambda kv: (-kv[1], kv[0]))
    (top_base, top_w), (_, second_w) = ranked[0], ranked[1]
    rest_w = sum(w for _, w in ranked[1:])

    if coverage < MIN_COVERAGE:
        return "N", top_w, coverage, weights, False
    # Floating point safe "equal weight" tie detection.
    tie = second_w > 0.0 and abs(top_w - second_w) <= 1e-9 * max(top_w, 1.0)
    if tie:
        return "N", top_w, coverage, weights, True
    if top_w < WEIGHT_MARGIN * rest_w:
        return "N", top_w, coverage, weights, False
    return top_base, top_w, coverage, weights, False


def build_consensus(
    req: ConsensusRequest,
) -> ConsensusResult:
    """Run the whole pipeline over a validated request (version assigned later)."""
    reads, duplicate_rows = dedupe_reads(req)

    by_locus: Dict[str, List[ReadIn]] = defaultdict(list)
    for r in reads:
        by_locus[r.locus].append(r)

    consensus_out: List[dict] = []
    clusters_out: List[dict] = []
    eliminated_out: List[dict] = []
    n_ambiguous = 0

    for locus in sorted(by_locus):
        locus_reads = by_locus[locus]
        reads_by_umi: Dict[str, List[ReadIn]] = defaultdict(list)
        for r in locus_reads:
            reads_by_umi[r.umi].append(r)

        for cl in cluster_locus(locus_reads):
            cl_reads: List[ReadIn] = []
            read_ids: List[str] = []
            for umi in cl.members:
                for r in reads_by_umi[umi]:
                    cl_reads.append(r)
                    read_ids.append(r.read_id)
            depth = len(cl_reads)

            clusters_out.append(
                {
                    "locus": locus,
                    "representative_umi": cl.representative,
                    "members": sorted(set(cl.members)),
                    "ambiguous": cl.ambiguous,
                    "depth": depth,
                    "read_ids": sorted(read_ids),
                }
            )
            if cl.ambiguous:
                n_ambiguous += 1
                eliminated_out.append(
                    {
                        "locus": locus,
                        "representative_umi": cl.representative,
                        "reason": "UMI hit more than one cluster representative "
                        "(Hamming distance <= 1); whole cluster is ambiguous",
                        "reason_code": "ambiguous_multi_representative",
                        "ambiguous": True,
                        "depth": depth,
                        "n_positions": [],
                        "details": {"representatives": sorted(set(cl.seeds))},
                    }
                )
                continue

            seq_len = max(len(r.sequence) for r in cl_reads)
            calls: List[str] = []
            n_positions: List[int] = []
            position_info: List[dict] = []
            for pos in range(seq_len):
                base, top_w, coverage, weights, tie = _vote_position(cl_reads, pos)
                calls.append(base)
                if base == "N":
                    n_positions.append(pos)
                    position_info.append(
                        {
                            "position": pos,
                            "coverage": coverage,
                            "tie": tie,
                            "weights": {
                                b: round(w, 4) for b, w in weights.items()
                            },
                        }
                    )

            if n_positions:
                eliminated_out.append(
                    {
                        "locus": locus,
                        "representative_umi": cl.representative,
                        "reason": f"{len(n_positions)} position(s) failed coverage "
                        "or Phred-weight dominance and were called N",
                        "reason_code": "consensus_contains_n",
                        "ambiguous": False,
                        "depth": depth,
                        "n_positions": n_positions,
                        "details": {"positions": position_info},
                    }
                )
            else:
                consensus_out.append(
                    {
                        "locus": locus,
                        "representative_umi": cl.representative,
                        "sequence": "".join(calls),
                        "depth": depth,
                        "ambiguous": False,
                    }
                )

    stats = {
        "loci": len(by_locus),
        "raw_reads": len(req.reads),
        "distinct_reads": len(reads),
        "duplicate_read_rows": duplicate_rows,
        "clusters": len(clusters_out),
        "consensus": len(consensus_out),
        "eliminated": len(eliminated_out),
        "ambiguous_clusters": n_ambiguous,
    }

    return {  # type: ignore[return-value]
        "version": 0,
        "request_id": req.request_id,
        "sample_id": req.sample_id,
        "consensus": consensus_out,
        "clusters": clusters_out,
        "eliminated": eliminated_out,
        "stats": stats,
    }
