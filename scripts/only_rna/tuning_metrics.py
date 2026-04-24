from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CandidateScore:
    qc_score: float
    annotation_score: float
    embedding_score: float
    total_score: float
    reason_code: str


def _clamp_score(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def score_qc_metrics(*, n_cells_total: int, n_cells_pass_qc: int) -> float:
    if int(n_cells_total) <= 0:
        return 0.0
    return _clamp_score(float(n_cells_pass_qc) / float(n_cells_total))


def score_annotation_metrics(
    *,
    method_status: str,
    confidence_mean: float,
    low_confidence_fraction: float,
) -> float:
    if method_status != "ok":
        return 0.0
    score = float(confidence_mean) * (1.0 - float(low_confidence_fraction))
    return _clamp_score(score)


def score_embedding_metrics(
    *,
    separation_score: float,
    fragmentation_penalty: float,
) -> float:
    score = float(separation_score) - float(fragmentation_penalty)
    return _clamp_score(score)


def summarize_candidate_score(
    *,
    qc_score: float,
    annotation_score: float,
    embedding_score: float,
    reason_code: str,
) -> CandidateScore:
    total_score = (
        float(qc_score) + float(annotation_score) + float(embedding_score)
    ) / 3.0
    return CandidateScore(
        qc_score=float(qc_score),
        annotation_score=float(annotation_score),
        embedding_score=float(embedding_score),
        total_score=total_score,
        reason_code=reason_code,
    )


__all__ = [
    "CandidateScore",
    "score_annotation_metrics",
    "score_embedding_metrics",
    "score_qc_metrics",
    "summarize_candidate_score",
]
