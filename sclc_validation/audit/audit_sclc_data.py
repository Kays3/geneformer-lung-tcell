#!/usr/bin/env python3
"""Reproducible feasibility audit for open SCLC single-cell and spatial data."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import pickle
from pathlib import Path

import pandas as pd
from PIL import Image


AUDIT_DIR = Path(__file__).resolve().parent
SOURCE_DIR = AUDIT_DIR / "source_metadata"
RESULTS_DIR = AUDIT_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SIGNATURE_GENES = [
    "PDCD1",
    "CTLA4",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "TOX",
    "LAYN",
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMH",
    "IFNG",
    "TCF7",
    "SLAMF6",
    "IL7R",
    "CCR7",
    "ASCL1",
    "NEUROD1",
    "POU2F3",
    "YAP1",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tissue_compartment(label: str) -> str:
    if label.startswith("N"):
        return "adjacent_normal"
    if label == "M":
        return "relapse_or_metastatic"
    return "primary_tumor"


def audit_omix() -> tuple[dict, pd.DataFrame]:
    metadata_path = SOURCE_DIR / "OMIX002441-02.csv"
    matrix_path = SOURCE_DIR / "OMIX002441-01.csv"
    token_path = SOURCE_DIR / "token_dictionary_gc95M.pkl"
    name_map_path = SOURCE_DIR / "gene_name_id_dict_gc95M.pkl"

    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata["compartment"] = metadata["tissue"].map(tissue_compartment)

    with matrix_path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        genes = [row[0] for row in reader]

    with token_path.open("rb") as handle:
        token_dictionary = pickle.load(handle)
    with name_map_path.open("rb") as handle:
        gene_name_to_id = pickle.load(handle)

    unique_genes = set(genes)
    mapped_genes = unique_genes.intersection(gene_name_to_id)
    tokenable_genes = {
        gene
        for gene in mapped_genes
        if gene_name_to_id[gene] in token_dictionary
    }

    t_cells = metadata.loc[metadata["cell_type"].eq("Tcell")].copy()
    donor_compartment = pd.crosstab(
        t_cells["patient"], t_cells["compartment"]
    ).reindex(
        columns=["primary_tumor", "adjacent_normal", "relapse_or_metastatic"],
        fill_value=0,
    )
    donor_compartment.index.name = "patient"
    donor_compartment.to_csv(RESULTS_DIR / "omix_tcell_counts_by_patient.csv")

    qc = (
        t_cells.groupby("compartment")
        .agg(
            cells=("patient", "size"),
            patients=("patient", "nunique"),
            median_features=("nFeature_RNA", "median"),
            median_counts=("nCount_RNA", "median"),
            median_mapping_rate=("MappingRate", "median"),
        )
        .reset_index()
    )
    qc.to_csv(RESULTS_DIR / "omix_tcell_qc_by_compartment.csv", index=False)

    paired_20 = int(
        (
            (donor_compartment["primary_tumor"] >= 20)
            & (donor_compartment["adjacent_normal"] >= 20)
        ).sum()
    )
    result = {
        "accession": "OMIX002441",
        "matrix_cells": len(header) - 1,
        "metadata_cells": int(len(metadata)),
        "matrix_metadata_exact_order_match": header[1:] == metadata.index.tolist(),
        "patients": int(metadata["patient"].nunique()),
        "t_cells": int(len(t_cells)),
        "primary_tumor_t_cells": int(
            t_cells["compartment"].eq("primary_tumor").sum()
        ),
        "adjacent_normal_t_cells": int(
            t_cells["compartment"].eq("adjacent_normal").sum()
        ),
        "relapse_or_metastatic_t_cells": int(
            t_cells["compartment"].eq("relapse_or_metastatic").sum()
        ),
        "patients_with_20_primary_t_cells": int(
            (donor_compartment["primary_tumor"] >= 20).sum()
        ),
        "patients_with_20_adjacent_normal_t_cells": int(
            (donor_compartment["adjacent_normal"] >= 20).sum()
        ),
        "patients_with_20_t_cells_in_both": paired_20,
        "metadata_null_cells": int(metadata.isna().sum().sum()),
        "duplicate_cell_ids": int(metadata.index.duplicated().sum()),
        "genes": len(genes),
        "unique_genes": len(unique_genes),
        "duplicate_gene_names": len(genes) - len(unique_genes),
        "geneformer_name_mapped_genes": len(mapped_genes),
        "geneformer_name_mapping_rate": len(mapped_genes) / len(unique_genes),
        "geneformer_tokenable_genes": len(tokenable_genes),
        "geneformer_tokenable_rate": len(tokenable_genes) / len(unique_genes),
        "signature_genes_present": {
            gene: gene in unique_genes for gene in SIGNATURE_GENES
        },
        "signature_genes_tokenable": {
            gene: (
                gene in gene_name_to_id
                and gene_name_to_id[gene] in token_dictionary
            )
            for gene in SIGNATURE_GENES
        },
        "matrix_sha256": sha256(matrix_path),
        "metadata_sha256": sha256(metadata_path),
        "token_dictionary_sha256": sha256(token_path),
        "gene_name_map_sha256": sha256(name_map_path),
    }
    return result, donor_compartment


def matrix_market_dimensions(path: Path) -> tuple[str, int, int, int]:
    with gzip.open(path, "rt") as handle:
        banner = handle.readline().strip()
        for line in handle:
            if not line.startswith("%"):
                rows, columns, nonzero = map(int, line.split())
                return banner, rows, columns, nonzero
    raise ValueError(f"No Matrix Market dimensions found in {path}")


def audit_spatial() -> tuple[list[dict], pd.DataFrame]:
    raw_dir = SOURCE_DIR / "GSE263196_RAW"
    rows: list[dict] = []
    for matrix_path in sorted(raw_dir.glob("*_matrix.mtx.gz")):
        prefix = matrix_path.name.removesuffix("_matrix.mtx.gz")
        barcode_path = raw_dir / f"{prefix}_barcodes.tsv.gz"
        feature_path = raw_dir / f"{prefix}_features.tsv.gz"
        position_path = raw_dir / f"{prefix}_tissue_positions_list.csv.gz"
        image_path = raw_dir / f"{prefix}_tissue_lowres_image.png.gz"
        scalefactor_path = raw_dir / f"{prefix}_scalefactors_json.json.gz"

        with gzip.open(barcode_path, "rt") as handle:
            barcodes = [line.rstrip("\n") for line in handle]
        with gzip.open(feature_path, "rt") as handle:
            features = [line.rstrip("\n").split("\t") for line in handle]
        with gzip.open(position_path, "rt") as handle:
            positions = list(csv.reader(handle))
        with gzip.open(image_path, "rb") as handle:
            with Image.open(io.BytesIO(handle.read())) as image:
                image_width, image_height = image.size
        with gzip.open(scalefactor_path, "rt") as handle:
            scalefactors = json.load(handle)

        banner, feature_rows, barcode_columns, nonzero = (
            matrix_market_dimensions(matrix_path)
        )
        feature_symbols = {row[1] for row in features}
        position_barcodes = {row[0] for row in positions}
        in_tissue_spots = sum(int(row[1]) for row in positions)
        rows.append(
            {
                "sample": prefix,
                "matrix_banner": banner,
                "features": len(features),
                "barcodes": len(barcodes),
                "positions": len(positions),
                "in_tissue_spots": in_tissue_spots,
                "matrix_feature_rows": feature_rows,
                "matrix_barcode_columns": barcode_columns,
                "matrix_nonzero_values": nonzero,
                "features_match_matrix": len(features) == feature_rows,
                "barcodes_match_matrix": len(barcodes) == barcode_columns,
                "positions_cover_matrix_barcodes": set(barcodes).issubset(
                    position_barcodes
                ),
                "unique_ensembl_ids": len({row[0] for row in features}),
                "unique_gene_symbols": len(feature_symbols),
                "signature_genes_present": sum(
                    gene in feature_symbols for gene in SIGNATURE_GENES
                ),
                "signature_genes_expected": len(SIGNATURE_GENES),
                "image_width": image_width,
                "image_height": image_height,
                "spot_diameter_fullres": scalefactors.get(
                    "spot_diameter_fullres"
                ),
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(RESULTS_DIR / "gse263196_spatial_file_audit.csv", index=False)
    return rows, table


def audit_geomx() -> tuple[dict, pd.DataFrame]:
    normalized_path = (
        SOURCE_DIR / "GSE261348_IMfirst_DSP_normalizedcounts.xlsx"
    )
    raw_path = SOURCE_DIR / "GSE261348_IMfirst_DSP_rawcounts.xlsx"

    segments = pd.read_csv(
        SOURCE_DIR / "GSE261348_normalizedcounts_SegmentProperties.csv"
    )
    counts = pd.read_csv(
        SOURCE_DIR / "GSE261348_normalizedcounts_TargetCountMatrix.csv"
    )
    targets = set(counts.iloc[:, 0].astype(str))
    matrix_segments = list(map(str, counts.columns[1:]))
    property_segments = set(segments["SegmentDisplayName"].astype(str))
    excluded_segments = sorted(property_segments.difference(matrix_segments))
    qc_flags = segments["QCFlags"].dropna().astype(str)

    raw_counts = pd.read_csv(
        SOURCE_DIR / "GSE261348_rawcounts_TargetCountMatrix.csv"
    )
    numeric = raw_counts.iloc[:, 1:].to_numpy(dtype=float)
    raw_matrix_is_nonnegative_integer = bool(
        (numeric >= 0).all() and (numeric == numeric.round()).all()
    )
    clinical_keywords = (
        "response",
        "outcome",
        "survival",
        "benefit",
        "progression",
        "death",
        "censor",
        "treatment",
    )
    outcome_columns = [
        column
        for column in map(str, segments.columns)
        if any(keyword in column.lower() for keyword in clinical_keywords)
    ]
    segment_export = segments[
        [
            "SegmentDisplayName",
            "SlideName",
            "ROILabel",
            "SegmentLabel",
            "QCFlags",
            "NormalizationFactor",
        ]
    ].copy()
    segment_export["included_in_expression_matrix"] = segment_export[
        "SegmentDisplayName"
    ].astype(str).isin(matrix_segments)
    segment_export.to_csv(
        RESULTS_DIR / "gse261348_segment_audit.csv", index=False
    )
    result = {
        "accession": "GSE261348",
        "segment_property_rows": int(len(segments)),
        "expression_aois": len(matrix_segments),
        "unique_slides": int(segments["SlideName"].nunique()),
        "targets": len(targets),
        "qc_flagged_segments": int(segments["QCFlags"].notna().sum()),
        "qc_flag_values": sorted(qc_flags.unique()),
        "segments_excluded_from_expression_matrix": len(excluded_segments),
        "excluded_segment_names": excluded_segments,
        "signature_genes_present": {
            gene: gene in targets for gene in SIGNATURE_GENES
        },
        "signature_genes_present_count": sum(
            gene in targets for gene in SIGNATURE_GENES
        ),
        "outcome_columns_in_workbook": outcome_columns,
        "raw_workbook_matrix_is_nonnegative_integer": (
            raw_matrix_is_nonnegative_integer
        ),
        "normalized_workbook_sha256": sha256(normalized_path),
        "raw_workbook_sha256": sha256(raw_path),
    }
    return result, segment_export


def audit_cellxgene_collection() -> dict:
    path = SOURCE_DIR / "cellxgene_collection_62e8f058.json"
    collection = json.loads(path.read_text())
    datasets = []
    for dataset in collection["datasets"]:
        asset = next(
            asset
            for asset in dataset["assets"]
            if asset["filetype"].upper() == "H5AD"
        )
        datasets.append(
            {
                "title": dataset["title"],
                "cell_count": dataset["cell_count"],
                "feature_count": dataset["feature_count"],
                "donor_count": len(dataset["donor_id"]),
                "diseases": sorted(item["label"] for item in dataset["disease"]),
                "raw_data_location": dataset.get("raw_data_location"),
                "h5ad_bytes": asset["filesize"],
                "h5ad_url": asset["url"],
                "dataset_id": dataset["dataset_id"],
            }
        )
    pd.DataFrame(datasets).to_csv(
        RESULTS_DIR / "cellxgene_collection_inventory.csv", index=False
    )
    t_cell = next(row for row in datasets if row["title"] == "T cells")
    return {
        "collection_id": collection["collection_id"],
        "description": collection["description"],
        "dataset_count": len(datasets),
        "t_cell_dataset": t_cell,
        "collection_metadata_sha256": sha256(path),
    }


def audit_htan_tcell_h5ad() -> tuple[dict | None, pd.DataFrame | None]:
    path = SOURCE_DIR / "HTAN_MSK_T_cells_6fde3ad9.h5ad"
    if not path.exists():
        return None, None

    import anndata as ad
    import numpy as np
    from scipy import sparse

    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs.copy()
    disease_column = next(
        column
        for column in ["disease", "disease_ontology_term_id"]
        if column in obs
    )
    donor_column = next(
        column
        for column in ["donor_id", "patient", "individual"]
        if column in obs
    )
    cell_type_column = next(
        column
        for column in ["cell_type", "cell_type_ontology_term_id"]
        if column in obs
    )

    disease_donor = (
        obs.groupby([disease_column, donor_column], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    disease_donor.to_csv(
        RESULTS_DIR / "htan_tcell_counts_by_disease_donor.csv", index=False
    )
    disease_summary = (
        disease_donor.groupby(disease_column, observed=True)
        .agg(
            cells=("cells", "sum"),
            donors=(donor_column, "nunique"),
            median_cells_per_donor=("cells", "median"),
            min_cells_per_donor=("cells", "min"),
            max_cells_per_donor=("cells", "max"),
            donors_with_20_cells=("cells", lambda values: int((values >= 20).sum())),
            donors_with_100_cells=(
                "cells",
                lambda values: int((values >= 100).sum()),
            ),
        )
        .reset_index()
    )
    disease_summary.to_csv(
        RESULTS_DIR / "htan_tcell_summary_by_disease.csv", index=False
    )

    raw_available = adata.raw is not None
    raw_integer_sample = None
    if raw_available:
        sample = adata.raw.X[: min(100, adata.n_obs), : min(1000, adata.raw.n_vars)]
        if sparse.issparse(sample):
            values = sample.data
        else:
            values = np.asarray(sample).ravel()
        raw_integer_sample = bool(
            np.all(values >= 0) and np.allclose(values, np.round(values))
        )

    var_names = set(map(str, adata.var_names))
    var_columns = list(map(str, adata.var.columns))
    raw_var = adata.raw.var if raw_available else None
    raw_var_columns = list(map(str, raw_var.columns)) if raw_available else []
    symbol_column = next(
        (
            column
            for column in [
                "feature_name",
                "gene_name",
                "gene_symbol",
                "symbol",
            ]
            if column in raw_var_columns
        ),
        None,
    )
    raw_symbols = (
        set(raw_var[symbol_column].astype(str))
        if symbol_column is not None
        else set()
    )
    result = {
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "obs_columns": list(map(str, obs.columns)),
        "disease_column": disease_column,
        "donor_column": donor_column,
        "cell_type_column": cell_type_column,
        "disease_values": sorted(map(str, obs[disease_column].unique())),
        "cell_type_values": sorted(map(str, obs[cell_type_column].unique())),
        "duplicate_cell_ids": int(obs.index.duplicated().sum()),
        "obs_null_cells": int(obs.isna().sum().sum()),
        "raw_available": raw_available,
        "raw_integer_sample": raw_integer_sample,
        "var_columns": var_columns,
        "raw_var_columns": raw_var_columns,
        "raw_gene_symbol_column": symbol_column,
        "var_names_are_ensembl_rate": sum(
            name.startswith("ENSG") for name in var_names
        )
        / len(var_names),
        "signature_genes_as_var_names": {
            gene: gene in var_names for gene in SIGNATURE_GENES
        },
        "signature_genes_in_raw_symbols": {
            gene: gene in raw_symbols for gene in SIGNATURE_GENES
        },
        "sha256": sha256(path),
    }
    adata.file.close()
    return result, disease_summary


def main() -> None:
    omix, _ = audit_omix()
    spatial_rows, _ = audit_spatial()
    geomx, _ = audit_geomx()
    collection = audit_cellxgene_collection()
    htan, _ = audit_htan_tcell_h5ad()

    summary = {
        "audit_version": "2026-07-28",
        "omix002441": omix,
        "gse263196": {
            "samples": len(spatial_rows),
            "all_matrix_integrity_checks_pass": all(
                row["features_match_matrix"]
                and row["barcodes_match_matrix"]
                and row["positions_cover_matrix_barcodes"]
                for row in spatial_rows
            ),
            "total_barcodes": sum(row["barcodes"] for row in spatial_rows),
            "total_in_tissue_spots": sum(
                row["in_tissue_spots"] for row in spatial_rows
            ),
            "all_signature_genes_present": all(
                row["signature_genes_present"]
                == row["signature_genes_expected"]
                for row in spatial_rows
            ),
            "archive_sha256": sha256(
                SOURCE_DIR / "GSE263196_RAW.tar"
            ),
        },
        "gse261348": geomx,
        "cellxgene_collection": collection,
        "htan_tcell_h5ad": htan,
    }
    (RESULTS_DIR / "audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
