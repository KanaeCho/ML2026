#!/usr/bin/env python3
"""Export co-ATAC sample matrix directories to sample-level h5ad files."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
import shutil


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_ROOT = ROOT / "output" / "co" / "atac"
DEFAULT_DATA_ROOT = ROOT / "data"


ATAC_OBS_COLUMNS = [
    "cell_barcode",
    "sample",
    "dataset",
    "age",
    "health",
    "donor",
    "final_celltype",
    "umap_atac_1",
    "umap_atac_2",
    "cima_ref_umap_1",
    "cima_ref_umap_2",
    "cima_cell_type_l1",
    "cima_cell_type_l1_masked",
    "cima_cell_type_l1_cluster_consensus",
    "cima_cell_type_l2",
    "cima_cell_type_l3",
    "cima_cell_type_l4",
    "cima_l1_low_confidence",
    "cima_l1_cluster_purity",
    "cima_l4_score",
    "cima_l4_score_margin",
]


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t")[0] for line in handle]


def sanitize_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.index = pd.Index(["" if pd.isna(v) else str(v) for v in out.index], dtype=object)
    out.columns = pd.Index([str(v) for v in out.columns], dtype=object)
    for column in out.columns:
        series = out[column]
        if (
            pd.api.types.is_string_dtype(series.dtype)
            or pd.api.types.is_object_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            out[column] = pd.Series(
                ["" if pd.isna(v) else str(v) for v in series.tolist()],
                index=series.index,
                dtype=object,
            )
    return out


def find_existing(path: Path, names: list[str]) -> Path:
    for name in names:
        candidate = path / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing one of {names} under {path}")


def co_metadata_lookup(data_root: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = data_root / "reference" / "co.xlsx"
    if not path.exists():
        return {}
    df = pd.read_excel(path, dtype=str).fillna("")
    required = {"sample", "dataset", "age", "health", "donor"}
    if not required <= set(df.columns):
        return {}
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for _, row in df.iterrows():
        dataset = str(row.get("dataset", "")).strip()
        sample = str(row.get("sample", "")).strip()
        if not dataset or not sample:
            continue
        lookup[(dataset, sample)] = {
            "sample": sample,
            "dataset": dataset,
            "age": str(row.get("age", "")).strip(),
            "health": str(row.get("health", "")).strip(),
            "donor": str(row.get("donor", "")).strip(),
        }
    return lookup


def export_sample(
    sample_dir: Path,
    *,
    overwrite: bool = False,
    data_root: Path = DEFAULT_DATA_ROOT,
    cleanup: bool = False,
) -> Path:
    sample_id = sample_dir.name
    output_path = sample_dir / f"{sample_id}.h5ad"
    if output_path.exists() and not overwrite:
        return output_path

    matrix_dir = sample_dir / "matrix"
    matrix_path = matrix_dir / "matrix.mtx"
    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)

    barcodes_path = find_existing(matrix_dir, ["barcodes.tsv.gz", "barcodes.tsv"])
    features_path = find_existing(matrix_dir, ["features.tsv.gz", "features.tsv"])
    metadata_path = sample_dir / "metadata.csv"
    metadata_qc_path = sample_dir / "metadata_qc.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)

    barcodes = read_lines(barcodes_path)
    features = read_lines(features_path)
    matrix = sparse.csr_matrix(mmread(matrix_path))
    if matrix.shape == (len(features), len(barcodes)):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (len(barcodes), len(features)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match {len(barcodes)} barcodes and {len(features)} features"
        )
    matrix = sparse.csr_matrix(matrix, dtype=np.float32)

    obs = pd.read_csv(metadata_path, dtype={"cell_barcode": str})
    if "cell_barcode" not in obs.columns:
        raise KeyError(f"{metadata_path} must contain cell_barcode")
    obs = obs.set_index("cell_barcode", drop=False)
    missing = [barcode for barcode in barcodes if barcode not in obs.index]
    if missing:
        raise ValueError(f"{sample_dir} metadata.csv is missing {len(missing)} matrix barcodes")
    obs = obs.loc[barcodes]

    if metadata_qc_path.exists():
        qc_obs = pd.read_csv(metadata_qc_path, dtype={"cell_barcode": str})
        qc_barcodes = qc_obs["cell_barcode"].astype(str).tolist()
    elif "pass_qc" in obs.columns:
        qc_barcodes = obs.index[obs["pass_qc"].fillna(False).astype(bool)].astype(str).tolist()
    else:
        qc_barcodes = barcodes

    qc_indices = [i for i, barcode in enumerate(barcodes) if barcode in set(qc_barcodes)]
    barcodes_qc = [barcodes[i] for i in qc_indices]
    matrix = matrix[qc_indices, :].tocsr()
    obs = obs.loc[barcodes_qc].copy()

    dataset = str(obs["gse"].iloc[0]) if "gse" in obs.columns and len(obs) else sample_dir.parent.name
    meta = co_metadata_lookup(data_root).get((dataset, sample_id), {})
    obs["sample"] = meta.get("sample", sample_id)
    obs["dataset"] = meta.get("dataset", dataset)
    obs["age"] = meta.get("age", "")
    obs["health"] = meta.get("health", "")
    obs["donor"] = meta.get("donor", str(obs.get("individual_id", pd.Series([""])).iloc[0]) if len(obs) else "")
    obs["final_celltype"] = obs["cima_cell_type_l1"].astype(str) if "cima_cell_type_l1" in obs.columns else ""
    for column in ATAC_OBS_COLUMNS:
        if column not in obs.columns:
            obs[column] = ""
    obs = obs[ATAC_OBS_COLUMNS].copy()

    var = pd.DataFrame(index=pd.Index(features, name="feature_id"))
    var["feature_id"] = features
    var["feature_name"] = features

    adata = ad.AnnData(X=matrix, obs=sanitize_dataframe(obs), var=sanitize_dataframe(var))
    adata.uns["source_output_dir"] = str(sample_dir)
    adata.uns["gse"] = str(obs["dataset"].iloc[0]) if "dataset" in obs.columns and len(obs) else dataset
    adata.uns["sample_id"] = sample_id
    adata.uns["modality"] = "ATAC"
    adata.write_h5ad(output_path, convert_strings_to_categoricals=False)
    if cleanup:
        if matrix_dir.exists():
            shutil.rmtree(matrix_dir)
        rds_path = sample_dir / f"{sample_id}_seurat_qc.rds"
        if rds_path.exists():
            rds_path.unlink()
    return output_path


def discover_samples(input_root: Path) -> list[Path]:
    return sorted(path for path in input_root.glob("*/*") if path.is_dir() and (path / "metadata.csv").exists())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--sample-dir", type=Path, default=None)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--cleanup", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    samples = [args.sample_dir] if args.sample_dir is not None else discover_samples(args.input_root)
    if not samples:
        raise SystemExit(f"No co-ATAC sample directories found under {args.input_root}")
    written = []
    for sample_dir in samples:
        written.append(
            export_sample(
                sample_dir,
                overwrite=args.overwrite,
                data_root=args.data_root,
                cleanup=args.cleanup,
            )
        )
    print(f"exported_h5ad={len(written)}")
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
