#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from typing import Any
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"

RANK_WEIGHTS = [1.0, 0.9, 0.8, 0.72, 0.64, 0.56, 0.5, 0.44, 0.39, 0.35]

BROAD_MARKERS = {
    "B cell lineage": {
        "MS4A1",
        "CD79A",
        "CD79B",
        "TCL1A",
        "BANK1",
        "PAX5",
        "BLK",
        "IGLL5",
        "CD22",
        "BACH2",
        "FCHO1",
        "APOBEC3G",
        "FAM46C",
        "MZB1",
        "JCHAIN",
        "TNFRSF17",
    },
    "T / NK": {
        "CD3D",
        "CD3E",
        "TRBC1",
        "TRBC2",
        "CD8A",
        "CD8B",
        "IL7R",
        "LEF1",
        "CCR7",
        "BCL11B",
        "CRTAM",
        "PRKCQ",
        "IL2RB",
        "TNFRSF4",
        "TNFSF9",
        "PRF1",
        "GNLY",
        "NKG7",
        "KLRD1",
        "KLRF1",
        "CTSW",
        "UNC13D",
        "KIR2DL4",
    },
    "Myeloid / innate-like": {
        "CD14",
        "LYZ",
        "S100A8",
        "S100A9",
        "FCGR3A",
        "MS4A7",
        "LST1",
        "CST3",
        "VCAN",
        "FCN1",
        "APOBEC3A",
        "RORC",
        "PRG4",
        "PFKFB3",
        "CALD1",
        "ADTRP",
        "SULF2",
        "IL1RL2",
        "CX3CR1",
        "F2R",
        "INPP5D",
        "FCER1G",
        "TYROBP",
    },
    "Dendritic cell": {
        "FCER1A",
        "HLA-DRA",
        "HLA-DQA1",
        "CLEC10A",
        "XCR1",
        "CLEC9A",
        "BATF3",
        "IRF8",
        "LILRA4",
        "CLEC4C",
        "CADM1",
    },
    "Megakaryocyte / platelet": {
        "PPBP",
        "PF4",
        "GP1BA",
        "GP9",
        "GP5",
        "ITGA2B",
        "TUBB1",
        "PF4V1",
    },
}

SUBTYPE_MARKERS = {
    "B_CD79A_BLK_like": {"CD79A", "CD79B", "BLK", "PAX5", "MS4A1", "BANK1", "TCL1A"},
    "B_IGLL5_CECR3_like": {"IGLL5", "CECR3", "FAM46C", "MZB1", "JCHAIN", "TNFRSF17"},
    "CD8_T_CRTAM_like": {"CD8A", "CD8B", "CRTAM", "CTSW", "CCL5", "GZMK"},
    "Cytotoxic_T_NK_like": {
        "PRF1",
        "GNLY",
        "NKG7",
        "KLRD1",
        "KIR2DL4",
        "UNC13D",
        "GZMB",
    },
    "T_NK_BCL11B_like": {
        "BCL11B",
        "CD3D",
        "CD3E",
        "TRBC1",
        "TRBC2",
    },
    "Activated_TNFSF9_like": {"TNFSF9", "CD8A", "PRKCQ", "IL7R"},
    "Activated_IL2RB_TNFRSF4_like": {"IL2RB", "TNFRSF4", "HDC"},
    "Monocyte_CALD1_PRG4_like": {"CALD1", "PRG4", "CD14", "VCAN", "FCN1", "LYZ"},
    "Innate_APOBEC3A_RORC_like": {"APOBEC3A", "RORC", "FCGR3A", "MS4A7", "LST1"},
    "Myeloid_PFKFB3_like": {"PFKFB3", "IL1B", "LYZ", "FCN1", "PRG4"},
    "Myeloid_ADTRP_SULF2_candidate": {"ADTRP", "SULF2", "IL1RL2"},
    "Myeloid_CX3CR1_F2R_like": {"CX3CR1", "F2R", "INPP5D", "FCER1G", "TYROBP"},
    "cDC1_like": {"XCR1", "CLEC9A", "BATF3", "IRF8"},
    "cDC2_like": {"FCER1A", "CLEC10A", "HLA-DRA", "HLA-DQA1", "CST3"},
    "Platelet_GP5_like": {"GP5", "GP1BA", "GP9", "PF4", "PPBP", "ITGA2B", "TUBB1"},
}

SUBTYPE_TO_BROAD = {
    "B_CD79A_BLK_like": "B cell lineage",
    "B_IGLL5_CECR3_like": "B cell lineage",
    "CD8_T_CRTAM_like": "T / NK",
    "Cytotoxic_T_NK_like": "T / NK",
    "T_NK_BCL11B_like": "T / NK",
    "Activated_TNFSF9_like": "T / NK",
    "Activated_IL2RB_TNFRSF4_like": "T / NK",
    "Monocyte_CALD1_PRG4_like": "Myeloid / innate-like",
    "Innate_APOBEC3A_RORC_like": "Myeloid / innate-like",
    "Myeloid_PFKFB3_like": "Myeloid / innate-like",
    "Myeloid_ADTRP_SULF2_candidate": "Myeloid / innate-like",
    "Myeloid_CX3CR1_F2R_like": "Myeloid / innate-like",
    "cDC1_like": "Dendritic cell",
    "cDC2_like": "Dendritic cell",
    "Platelet_GP5_like": "Megakaryocyte / platelet",
}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_cluster(cluster: str) -> str:
    return cluster[1:] if cluster.startswith("g") else cluster


def slugify_label(value: str) -> str:
    return (
        value.replace(" / ", "_")
        .replace("/", "_")
        .replace(" ", "_")
        .replace("-", "_")
        .lower()
    )


def weighted_score(genes: list[str], marker_set: set[str]) -> tuple[float, list[str]]:
    seen = []
    score = 0.0
    for idx, gene in enumerate(genes):
        if gene in marker_set:
            score += RANK_WEIGHTS[min(idx, len(RANK_WEIGHTS) - 1)]
            if gene not in seen:
                seen.append(gene)
    return score, seen


def top_two(score_map: dict[str, float]) -> tuple[tuple[str, float], tuple[str, float]]:
    ordered = sorted(score_map.items(), key=lambda item: (-item[1], item[0]))
    first = ordered[0] if ordered else ("", 0.0)
    second = ordered[1] if len(ordered) > 1 else ("", 0.0)
    return first, second


def euclidean(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


def broad_confidence(
    top_score: float, second_score: float, hit_count: int, mode: str
) -> str:
    if mode == "relaxed":
        if hit_count >= 2 and top_score >= 1.0 and top_score >= second_score + 0.35:
            return "high"
        if hit_count >= 1 and top_score >= 0.45 and top_score >= second_score + 0.05:
            return "medium"
        if hit_count >= 1 and top_score >= 0.35:
            return "low"
        return "unresolved"

    if hit_count >= 2 and top_score >= 1.2 and top_score >= second_score + 0.45:
        return "high"
    if hit_count >= 1 and top_score >= 0.7 and top_score >= second_score + 0.2:
        return "medium"
    if hit_count >= 1 and top_score >= 0.5 and top_score >= second_score + 0.1:
        return "low"
    return "unresolved"


def build_cluster_gene_lists(
    marker_rows: list[dict[str, str]], top_n: int
) -> dict[str, list[str]]:
    genes: dict[str, list[str]] = defaultdict(list)
    for row in marker_rows:
        cluster = normalize_cluster(row["cluster"])
        gene = (row.get("gene_name") or "").strip().upper()
        if not gene:
            continue
        if gene not in genes[cluster]:
            genes[cluster].append(gene)
    return {cluster: gene_list[:top_n] for cluster, gene_list in genes.items()}


def build_centroids(
    metadata_rows: list[dict[str, str]],
) -> dict[str, tuple[float, float]]:
    sums: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    for row in metadata_rows:
        cluster = row["seurat_clusters"]
        sums[cluster][0] += float(row["UMAP_Harmony_1"])
        sums[cluster][1] += float(row["UMAP_Harmony_2"])
        sums[cluster][2] += 1.0
    centroids = {}
    for cluster, (sx, sy, n) in sums.items():
        centroids[cluster] = (sx / n, sy / n)
    return centroids


def cluster_sizes(metadata_rows: list[dict[str, str]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in metadata_rows:
        counts[row["seurat_clusters"]] += 1
    return dict(counts)


def sample_dominance(sample_rows: list[dict[str, str]]) -> dict[str, tuple[str, float]]:
    out = {}
    if not sample_rows:
        return out
    header = [key for key in sample_rows[0].keys() if key != ""]
    for row in sample_rows:
        cluster = row[""]
        counts = {sample: int(row[sample]) for sample in header}
        total = sum(counts.values())
        if total == 0:
            out[cluster] = ("", 0.0)
            continue
        top_sample, top_count = max(counts.items(), key=lambda item: item[1])
        out[cluster] = (top_sample, top_count / total)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="De novo manual annotation from current cluster marker peaks."
    )
    parser.add_argument(
        "--analysis-dir",
        default=str(OUTPUT_ROOT / "integration_merged_without_unknown_analysis"),
    )
    parser.add_argument("--top-n-genes", type=int, default=10)
    parser.add_argument(
        "--mode", choices=["conservative", "relaxed"], default="conservative"
    )
    parser.add_argument("--output-suffix", default="_de_novo")
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    marker_path = analysis_dir / "cluster_top_accessible_peaks.csv"
    metadata_path = analysis_dir / "integrated_metadata.csv.gz"
    sample_count_path = analysis_dir / "cluster_by_sample_counts.csv"
    if not all(
        path.exists() for path in [marker_path, metadata_path, sample_count_path]
    ):
        raise SystemExit("Missing required analysis files for de novo annotation")

    marker_rows = read_csv_rows(marker_path)
    metadata_rows = read_csv_rows(metadata_path)
    sample_rows = read_csv_rows(sample_count_path)
    cluster_genes = build_cluster_gene_lists(marker_rows, args.top_n_genes)
    centroids = build_centroids(metadata_rows)
    sizes = cluster_sizes(metadata_rows)
    top_sample_fraction = sample_dominance(sample_rows)

    score_rows: list[dict[str, Any]] = []
    broad_assignments: dict[str, dict[str, Any]] = {}
    anchors: dict[str, str] = {}

    for cluster, genes in sorted(cluster_genes.items(), key=lambda item: int(item[0])):
        broad_scores = {}
        broad_hits = {}
        for broad, markers in BROAD_MARKERS.items():
            score, hits = weighted_score(genes, markers)
            broad_scores[broad] = score
            broad_hits[broad] = hits

        (top_broad, top_score), (_, second_score) = top_two(broad_scores)
        hit_count = len(broad_hits.get(top_broad, []))
        conf = broad_confidence(top_score, second_score, hit_count, args.mode)
        broad_assignments[cluster] = {
            "broad_label": top_broad
            if conf != "unresolved"
            else "Unknown / unresolved",
            "broad_confidence": conf if conf != "unresolved" else "low",
            "broad_score": round(top_score, 4),
            "broad_score_second": round(second_score, 4),
            "broad_hits": broad_hits.get(top_broad, []),
            "broad_hit_count": hit_count,
            "genes": genes,
        }
        if conf == "high" or (args.mode == "relaxed" and conf == "medium"):
            anchors[cluster] = top_broad

    for cluster, data in broad_assignments.items():
        if data["broad_label"] != "Unknown / unresolved":
            continue
        centroid = centroids[cluster]
        neighbor_votes: dict[str, float] = defaultdict(float)
        distances = []
        for other_cluster, broad in anchors.items():
            dist = euclidean(centroid, centroids[other_cluster])
            distances.append((dist, other_cluster, broad))
        distances.sort(key=lambda item: item[0])
        for dist, _, broad in distances[:3]:
            neighbor_votes[broad] += 1 / max(dist, 1e-6)
        allow_neighbor_rescue = (
            int(data.get("broad_hit_count", 0)) >= 1 or args.mode == "relaxed"
        )
        if neighbor_votes and allow_neighbor_rescue:
            (neighbor_label, neighbor_score), (_, neighbor_second) = top_two(
                dict(neighbor_votes)
            )
            dominance_ratio = 1.25 if args.mode == "relaxed" else 1.4
            if neighbor_score >= max(neighbor_second, 1e-6) * dominance_ratio:
                data["broad_label"] = neighbor_label
                data["broad_confidence"] = "low"
                data["neighbor_support"] = [
                    cluster_id
                    for _, cluster_id, broad in distances[:3]
                    if broad == neighbor_label
                ]
            else:
                data["neighbor_support"] = [
                    cluster_id for _, cluster_id, _ in distances[:3]
                ]
        else:
            data["neighbor_support"] = []

    annotation_rows: list[dict[str, object]] = []
    for cluster, data in sorted(
        broad_assignments.items(), key=lambda item: int(item[0])
    ):
        broad = data["broad_label"]
        broad_conf = data["broad_confidence"]
        genes = data["genes"]

        subtype_scores = {}
        subtype_hits = {}
        for subtype, markers in SUBTYPE_MARKERS.items():
            if SUBTYPE_TO_BROAD[subtype] != broad:
                continue
            score, hits = weighted_score(genes, markers)
            subtype_scores[subtype] = score
            subtype_hits[subtype] = hits

        if subtype_scores:
            (top_subtype, top_sub_score), (_, second_sub_score) = top_two(
                subtype_scores
            )
        else:
            top_subtype, top_sub_score, second_sub_score = "", 0.0, 0.0

        if broad == "Unknown / unresolved":
            subtype_label = "Unknown_unresolved"
            subtype_conf = "low"
            matched = []
        elif (
            top_sub_score >= (0.6 if args.mode == "relaxed" else 0.8)
            and top_sub_score
            >= second_sub_score + (0.1 if args.mode == "relaxed" else 0.2)
            and len(subtype_hits.get(top_subtype, [])) >= 1
        ):
            subtype_label = top_subtype
            subtype_conf = "high" if broad_conf == "high" else "medium"
            matched = subtype_hits[top_subtype]
        elif (
            top_sub_score >= (0.35 if args.mode == "relaxed" else 0.5)
            and top_sub_score
            >= second_sub_score + (0.05 if args.mode == "relaxed" else 0.15)
            and len(subtype_hits.get(top_subtype, [])) >= 1
        ):
            subtype_label = (
                top_subtype
                if top_subtype.endswith("_candidate")
                else f"{top_subtype}_candidate"
            )
            subtype_conf = "low"
            matched = subtype_hits[top_subtype]
        else:
            subtype_label = f"Unresolved_{slugify_label(broad)}"
            subtype_conf = "low"
            matched = data.get("broad_hits", [])

        if (
            broad == "Dendritic cell"
            and len(data.get("broad_hits", [])) < 2
            and top_sub_score < 0.5
        ):
            broad = "Unknown / unresolved"
            broad_conf = "low"
            subtype_label = "Unknown_unresolved"
            subtype_conf = "low"
            matched = data.get("broad_hits", [])

        top_sample, sample_fraction = top_sample_fraction.get(cluster, ("", 0.0))
        if sample_fraction >= 0.85 and sizes.get(cluster, 0) < 500:
            broad = "Unknown / unresolved"
            subtype_label = "Unknown_sample_biased"
            broad_conf = "low"
            subtype_conf = "low"

        evidence = matched[:3] if matched else genes[:3]
        annotation_confidence = (
            "high"
            if broad_conf == "high" and subtype_conf == "high"
            else (
                "medium"
                if broad_conf in {"high", "medium"}
                and subtype_conf in {"high", "medium"}
                else "low"
            )
        )

        annotation_rows.append(
            {
                "seurat_clusters": cluster,
                "celltype": broad,
                "celltype_subtype": subtype_label,
                "evidence_markers": "; ".join(evidence),
                "annotation_confidence": annotation_confidence,
            }
        )

        score_rows.append(
            {
                "seurat_clusters": cluster,
                "top_genes": "; ".join(genes),
                "broad_label": broad,
                "broad_confidence": broad_conf,
                "broad_score": data["broad_score"],
                "broad_score_second": data["broad_score_second"],
                "broad_hits": "; ".join(data.get("broad_hits", [])),
                "subtype_label": subtype_label,
                "subtype_score": round(top_sub_score, 4),
                "top_sample": top_sample,
                "top_sample_fraction": round(sample_fraction, 4),
                "neighbor_support": "; ".join(data.get("neighbor_support", [])),
            }
        )

    suffix = args.output_suffix
    map_path = analysis_dir / f"cluster_celltype_annotation_map{suffix}.csv"
    score_path = analysis_dir / f"cluster_marker_annotation_scores{suffix}.csv"
    notes_path = analysis_dir / f"celltype_annotation_notes{suffix}.md"

    write_csv(
        map_path,
        annotation_rows,
        [
            "seurat_clusters",
            "celltype",
            "celltype_subtype",
            "evidence_markers",
            "annotation_confidence",
        ],
    )
    write_csv(
        score_path,
        score_rows,
        [
            "seurat_clusters",
            "top_genes",
            "broad_label",
            "broad_confidence",
            "broad_score",
            "broad_score_second",
            "broad_hits",
            "subtype_label",
            "subtype_score",
            "top_sample",
            "top_sample_fraction",
            "neighbor_support",
        ],
    )

    note_lines = [
        "# De novo celltype annotation notes",
        "",
        "- This annotation does not use any previous manual labels or label transfer.",
        f"- Mode: `{args.mode}`.",
        "- Broad labels are assigned from current cluster marker peaks using canonical PBMC marker genes.",
        "- Subtypes are assigned only when current marker support is stronger; otherwise unresolved subtype labels are kept.",
        "- Small clusters with extreme single-sample dominance are forced to unresolved labels.",
        f"- `{score_path.name}` contains the raw marker-support summary used for this de novo call.",
    ]
    notes_path.write_text("\n".join(note_lines) + "\n", encoding="utf-8")

    print(f"[de-novo] wrote {map_path}")
    print(f"[de-novo] wrote {score_path}")
    print(f"[de-novo] wrote {notes_path}")


if __name__ == "__main__":
    main()
