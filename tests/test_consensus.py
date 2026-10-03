"""Unit tests for clustering and weighted voting (no HTTP layer)."""
import math

import pytest
from pydantic import ValidationError

from app.consensus import (
    build_consensus,
    cluster_locus,
    find_readid_conflicts,
    hamming_le1,
    payload_hash,
)
from app.schemas import ConsensusRequest, ReadIn


def r(locus, read_id, umi, seq, quals=None):
    if quals is None:
        quals = [30] * len(seq)
    return ReadIn(locus=locus, read_id=read_id, umi=umi, sequence=seq, qualities=quals)


def req(request_id, reads, sample_id="s1"):
    return ConsensusRequest(
        request_id=request_id, sample_id=sample_id, reads=list(reads)
    )


# --------------------------------------------------------------------------
# schema validation
# --------------------------------------------------------------------------

def test_qualities_must_match_sequence_length():
    with pytest.raises(ValidationError) as ei:
        r("L1", "x", "AAAA", "ACGT", [30, 30, 30])
    assert "qualities length" in str(ei.value)


@pytest.mark.parametrize("bad_q", [-1, 41])
def test_quality_range_enforced(bad_q):
    with pytest.raises(ValidationError):
        r("L1", "x", "AAAA", "ACGT", [30, bad_q, 30, 30])


def test_sequence_only_acgt():
    with pytest.raises(ValidationError):
        r("L1", "x", "AAAA", "ACGN")


def test_umi_only_acgt():
    with pytest.raises(ValidationError):
        r("L1", "x", "AAA2", "ACGT")


def test_quality_ascii_phred33_accepted():
    read = ReadIn(
        locus="L1", read_id="x", umi="AAAA", sequence="ACGT", qualities="????"
    )
    assert read.qualities == [30, 30, 30, 30]
    read2 = ReadIn(
        locus="L1", read_id="y", umi="AAAA", sequence="AC", qualities="?+"
    )
    assert read2.qualities == [30, 10]


def test_quality_digit_string_accepted():
    read = ReadIn(
        locus="L1", read_id="x", umi="AAAA", sequence="ACG", qualities="30, 20, 10"
    )
    assert read.qualities == [30, 20, 10]


def test_extra_fields_rejected():
    with pytest.raises(ValidationError):
        ConsensusRequest.model_validate(
            {
                "request_id": "q",
                "sample_id": "s",
                "reads": [
                    {
                        "locus": "L1",
                        "read_id": "a",
                        "umi": "AAAA",
                        "sequence": "ACGT",
                        "qualities": [30] * 4,
                        "extra": 1,
                    }
                ],
            }
        )


# --------------------------------------------------------------------------
# hamming / clustering
# --------------------------------------------------------------------------

def test_hamming_helper():
    assert hamming_le1("AAAA", "AAAA")
    assert hamming_le1("AAAA", "AAGA")
    assert not hamming_le1("AAAA", "AGGA")
    assert not hamming_le1("AAA", "AAAA")  # different length never matches


def _reps(clusters):
    return [(c.representative, sorted(c.members), c.ambiguous) for c in clusters]


def test_frequency_desc_lexicographic_chooses_representative():
    reads = []
    reads += [r("L", f"g{i}", "AAAAA", "ACGT") for i in range(3)]
    reads += [r("L", f"t{i}", "AAAAT", "ACGT") for i in range(5)]
    reads += [r("L", "c", "AAAAC", "ACGT")]
    clusters = cluster_locus(reads)
    assert len(clusters) == 1
    # AAAAT is most frequent -> representative, even though AAAAA is lex smaller
    assert clusters[0].representative == "AAAAT"
    assert clusters[0].ambiguous is False


def test_seeds_must_be_distance_2_plus():
    # AAAAA and AAAGG differ at 2 positions -> two seeds.
    reads = [r("L", "a", "AAAAA", "ACGT"), r("L", "b", "AAAGG", "ACGT")]
    clusters = cluster_locus(reads)
    assert len(clusters) == 2


def test_multi_representative_hit_merges_and_marks_ambiguous():
    reads = [r("L", "a", "AAAAA", "ACGT"), r("L", "b", "AAAGG", "ACGT")]
    # AAAAG is distance 1 from BOTH AAAAA and AAAAGG.
    reads.append(r("L", "bridge", "AAAAG", "ACGT"))
    clusters = cluster_locus(reads)
    assert len(clusters) == 1
    assert clusters[0].ambiguous is True
    assert sorted(clusters[0].seeds) == ["AAAAA", "AAAGG"]


def test_two_bridges_between_same_reps_do_not_duplicate_seeds():
    # Two representatives at Hamming distance 2:
    #   A=GGGGGG, B=GGGGTT (differ at positions 4 and 5)
    # Two distinct UMIs are distance 1 from BOTH representatives:
    #   GGGGGT (d1 from A at pos5; d1 from B at pos4)
    #   GGGGTG (d1 from A at pos4; d1 from B at pos5)
    reads = [
        r("L", "a1", "GGGGGG", "ACGT"),
        r("L", "a2", "GGGGGG", "ACGT"),
        r("L", "b1", "GGGGTT", "ACGT"),
        r("L", "b2", "GGGGTT", "ACGT"),
        r("L", "br1", "GGGGGT", "ACGT"),
        r("L", "br2", "GGGGTG", "ACGT"),
    ]
    clusters = cluster_locus(reads)
    assert len(clusters) == 1
    cl = clusters[0]
    assert cl.ambiguous is True
    # Seeds must be reported exactly once even though two UMIs bridge them.
    assert sorted(cl.seeds) == ["GGGGGG", "GGGGTT"]
    assert set(cl.members) == {"GGGGGG", "GGGGGT", "GGGGTT", "GGGGTG"}


def test_read_id_conflict_detected():
    payload = req(
        "q",
        [
            r("L", "same-id", "AAAAA", "ACGT"),
            r("L", "same-id", "AAAAT", "ACGT"),
        ],
    )
    conflicts = find_readid_conflicts(payload)
    assert len(conflicts) == 1
    assert conflicts[0]["read_id"] == "same-id"


# --------------------------------------------------------------------------
# weighted voting
# --------------------------------------------------------------------------

def test_unanimous_consensus_is_called():
    reads = [r("L", f"r{i}", "AAAAA", "ACGT") for i in range(4)]
    out = build_consensus(req("q", reads))
    assert out["consensus"][0]["sequence"] == "ACGT"
    assert out["stats"]["consensus"] == 1
    assert out["stats"]["eliminated"] == 0


def test_coverage_below_three_eliminates():
    reads = [r("L", f"r{i}", "AAAAA", "ACGT") for i in range(2)]
    out = build_consensus(req("q", reads))
    assert out["consensus"] == []
    el = out["eliminated"][0]
    assert el["reason_code"] == "consensus_contains_n"
    assert el["n_positions"] == [0, 1, 2, 3]


def test_weight_dominance_two_high_quality_beat_weak_third():
    # Three reads at position 0: 2x A@Q40 vs 1x C@Q0 => A weight 20000 vs 1.
    reads = [
        r("L", "a", "AAAAA", "ACGT", [40, 30, 30, 30]),
        r("L", "b", "AAAAA", "ACGT", [40, 30, 30, 30]),
        r("L", "c", "AAAAA", "CCGT", [0, 30, 30, 30]),
    ]
    out = build_consensus(req("q", reads))
    assert out["consensus"][0]["sequence"][0] == "A"


def test_low_quality_outvoted_when_margin_not_double():
    # 2x A@Q30 (2000) vs 1x C@30 (1000): 2000 < 2 * 1000 strictly false?
    # 2000 == 2*1000 -> condition top >= 2*rest passes (inclusive).
    reads = [
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "b", "AAAAA", "ACGT"),
        r("L", "c", "AAAAA", "CCGT"),
    ]
    out = build_consensus(req("q", reads))
    assert out["consensus"][0]["sequence"][0] == "A"

    # Now 2x A@Q27 (~501.2 each = 1002) vs 1x C@30 (1000):
    # 1002 < 2*1000 -> N
    reads2 = [
        r("L", "a", "AAAAA", "ACGT", [27, 30, 30, 30]),
        r("L", "b", "AAAAA", "ACGT", [27, 30, 30, 30]),
        r("L", "c", "AAAAA", "CCGT", [30, 30, 30, 30]),
    ]
    out2 = build_consensus(req("q2", reads2))
    assert out2["consensus"] == []
    assert out2["eliminated"][0]["n_positions"] == [0]


def test_exact_weight_tie_is_n():
    reads = [
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "b", "AAAAA", "ACGT"),
        r("L", "c", "AAAAA", "CCGT"),
        r("L", "d", "AAAAA", "CCGT"),
    ]
    out = build_consensus(req("q", reads))
    el = out["eliminated"][0]
    assert el["n_positions"] == [0]
    assert el["details"]["positions"][0]["tie"] is True


def test_ambiguous_cluster_is_eliminated_without_voting():
    reads = [
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "b", "AAAGG", "ACGT"),
        r("L", "bridge", "AAAAG", "ACGT"),
        r("L", "a2", "AAAAA", "ACGT"),
        r("L", "b2", "AAAGG", "ACGT"),
        r("L", "br2", "AAAAG", "ACGT"),
    ]
    out = build_consensus(req("q", reads))
    assert out["consensus"] == []
    el = out["eliminated"][0]
    assert el["reason_code"] == "ambiguous_multi_representative"
    assert el["depth"] == 6
    assert out["stats"]["ambiguous_clusters"] == 1


def test_duplicate_rows_deduped_and_counted():
    reads = [
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "b", "AAAAA", "ACGT"),
        r("L", "c", "AAAAA", "ACGT"),
    ]
    out = build_consensus(req("q", reads))
    assert out["stats"]["raw_reads"] == 4
    assert out["stats"]["distinct_reads"] == 3
    assert out["stats"]["duplicate_read_rows"] == 1
    assert out["consensus"][0]["depth"] == 3


def test_results_sorted_by_locus():
    reads = [r("Z", f"z{i}", "AAAAA", "ACGT") for i in range(3)]
    reads += [r("A", f"a{i}", "AAAAA", "ACGT") for i in range(3)]
    reads += [r("M", f"m{i}", "AAAAA", "ACGT") for i in range(3)]
    out = build_consensus(req("q", reads))
    loci = [c["locus"] for c in out["consensus"]]
    assert loci == ["A", "M", "Z"]


def test_payload_hash_order_independent():
    a = [
        r("L", "a", "AAAAA", "ACGT"),
        r("L", "b", "AAAAT", "ACGT"),
    ]
    b = list(reversed(a))
    assert payload_hash(req("q", a)) == payload_hash(req("q", b))


def test_payload_hash_detects_content_change():
    a = [r("L", "a", "AAAAA", "ACGT")]
    b = [r("L", "a", "AAAAA", "ACGN", None)] if False else [
        r("L", "a", "AAAAA", "ACGG")
    ]
    assert payload_hash(req("q", a)) != payload_hash(req("q", b))
