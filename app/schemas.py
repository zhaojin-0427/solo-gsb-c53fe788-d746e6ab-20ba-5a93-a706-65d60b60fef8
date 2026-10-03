"""Request / response schemas for the UMI consensus denoising API."""
from __future__ import annotations

from typing import List, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

_BASES = "ACGT"


class ReadIn(BaseModel):
    """One sequenced read carrying a UMI."""

    # Extra fields are rejected so a typo cannot silently change a payload hash.
    model_config = ConfigDict(extra="forbid")

    locus: str = Field(min_length=1)
    read_id: str = Field(min_length=1)
    umi: str = Field(min_length=1)
    sequence: str = Field(min_length=1)
    # Either a list of integer Q-values or an already-encoded quality string
    # (comma/space separated digits, or FASTQ ASCII Phred+33).
    qualities: Union[List[int], str]

    @field_validator("locus", "read_id", "umi", "sequence")
    @classmethod
    def _no_whitespace(cls, v: str) -> str:
        v = v.strip()
        if any(ch.isspace() for ch in v):
            raise ValueError("must not contain whitespace")
        return v

    @field_validator("umi")
    @classmethod
    def _umi_bases(cls, v: str) -> str:
        v = v.upper()
        bad = {ch for ch in v if ch not in _BASES}
        if bad:
            raise ValueError(f"UMI may only contain A/C/G/T, found {sorted(bad)}")
        return v

    @field_validator("sequence")
    @classmethod
    def _seq_bases(cls, v: str) -> str:
        v = v.upper()
        bad = {ch for ch in v if ch not in _BASES}
        if bad:
            raise ValueError(f"sequence may only contain A/C/G/T, found {sorted(bad)}")
        return v

    @field_validator("qualities")
    @classmethod
    def _normalize_qualities(cls, v):
        if isinstance(v, list):
            qs: List[int] = []
            for q in v:
                if isinstance(q, bool) or not isinstance(q, int):
                    raise ValueError("each quality value must be an integer")
                qs.append(q)
            return qs
        raw = v.strip()
        # Digits (optionally comma/space separated) are literal Q-values;
        # any non-digit printable character is treated as Phred+33 ASCII.
        tokens = [t for t in raw.replace(",", " ").split() if t != ""]
        if tokens and all(t.isdigit() for t in tokens):
            return [int(t) for t in tokens]
        return [ord(ch) - 33 for ch in raw]

    @model_validator(mode="after")
    def _check_ranges_and_length(self) -> "ReadIn":
        if len(self.sequence) != len(self.qualities):
            raise ValueError(
                f"sequence length {len(self.sequence)} != qualities length "
                f"{len(self.qualities)} for read_id={self.read_id!r}"
            )
        for i, q in enumerate(self.qualities):
            if not 0 <= q <= 40:
                raise ValueError(
                    f"quality value {q} at position {i} of read_id={self.read_id!r} "
                    "is outside the allowed range 0..40"
                )
        return self


class ConsensusRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    request_id: str = Field(min_length=1)
    sample_id: str = Field(min_length=1)
    reads: List[ReadIn] = Field(min_length=1)

    @field_validator("request_id", "sample_id")
    @classmethod
    def _strip(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("must not be empty")
        return v


class ConsensusOut(BaseModel):
    locus: str
    representative_umi: str
    sequence: str
    depth: int
    ambiguous: bool


class ClusterOut(BaseModel):
    locus: str
    representative_umi: str
    members: List[str]
    ambiguous: bool
    depth: int
    read_ids: List[str]


class EliminatedOut(BaseModel):
    locus: str
    representative_umi: str
    reason: str
    reason_code: Literal[
        "ambiguous_multi_representative", "consensus_contains_n"
    ]
    ambiguous: bool
    depth: int
    n_positions: List[int]
    details: dict


class StatsOut(BaseModel):
    loci: int
    raw_reads: int
    distinct_reads: int
    duplicate_read_rows: int
    clusters: int
    consensus: int
    eliminated: int
    ambiguous_clusters: int


class ConsensusResult(BaseModel):
    version: int
    request_id: str
    sample_id: str
    consensus: List[ConsensusOut]
    clusters: List[ClusterOut]
    eliminated: List[EliminatedOut]
    stats: StatsOut
