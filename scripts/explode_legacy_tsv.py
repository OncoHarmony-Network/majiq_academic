"""
explode_legacy_tsv.py

Post-process MAJIQ v2 tsvs to access relevant columns for comparison to new
quantifications

Author: Joseph K Aicher
"""

from typing import Dict

import pandas as pd


def _df_explode(
    df: pd.DataFrame, cols: Dict[str, type], sep: str = ";"
) -> pd.DataFrame:
    """Explode sep-delimited string columns of df and convert to specified type

    Parameters
    ----------
    df: pd.DataFrame
    cols: Dict[str, type]
        {col_name: col_type for each col_name that is sep-delimited to explode}
    sep: str
        Delimiter between entries in the columns
    """
    if not cols:
        return df
    exploded = (
        pd.DataFrame({x: df[x].str.split(sep) for x in cols.keys()})
        .explode(list(cols.keys()))
        .astype(cols)
    )
    return exploded.reset_index().join(df.drop(columns=list(cols.keys())), on="index")[
        df.columns
    ]


def _extract_coordinates_series(x: pd.Series) -> pd.DataFrame:
    """Assuming x is string encoding start-end"""
    return x.str.extract(r"(?P<start>\d+)-(?P<end>\d+)").astype(int)


def _extract_coordinates(df: pd.DataFrame, cols: Dict[str, str]) -> pd.DataFrame:
    """Extract coordinates, cols maps original column names to prefixes for start/end"""
    if not cols:
        return df
    concat_df = [
        df[x_orig].pipe(_extract_coordinates_series).add_prefix(x_prefix)
        for x_orig, x_prefix in cols.items()
    ]
    # add dataframe with non-coordinate columns
    concat_df.append(df.drop(columns=list(cols.keys())))
    # get new columns in order
    new_columns = [
        x
        for s in (
            (f"{cols[y]}start", f"{cols[y]}end") if y in cols else (y,)
            for y in df.columns
        )
        for x in s
    ]
    return pd.concat(concat_df, axis=1, join="inner")[new_columns]


def _extract_lsv_id_series(x: pd.Series) -> pd.DataFrame:
    df = x.str.extract(
        r"(?P<gene_id>.+):(?P<event_type>[st]):(?P<start>[^:-]+)-(?P<end>[^:-]+)"
    )
    return df.assign(
        start=lambda df: df.start.where(df.start != "na", "-1").astype(int),
        end=lambda df: df.end.where(df.end != "na", "-1").astype(int),
    ).rename(columns={x: f"ref_exon_{x}" for x in ("start", "end")})


def _extract_lsv_id(df: pd.DataFrame, lsv_id_col: str = "lsv_id") -> pd.DataFrame:
    df_lsv_id = _extract_lsv_id_series(df[lsv_id_col])
    # check overlapping columns
    for x in df.columns.intersection(df_lsv_id.columns):
        pd.testing.assert_series_equal(df[x], df_lsv_id[x])
    unique_lsv_columns = df_lsv_id.columns.difference(df.columns)
    # new columns in order
    new_columns = [
        x
        for s in (unique_lsv_columns if y == lsv_id_col else (y,) for y in df.columns)
        for x in s
    ]
    return pd.concat(
        [df_lsv_id[unique_lsv_columns], df.drop(columns=[lsv_id_col])],
        axis=1,
        join="inner",
    )[new_columns]


def extract_psi(df: pd.DataFrame) -> pd.DataFrame:
    """Extract columns from majiq (v2) psi.tsv"""
    PSI_EXPLODE_COLUMNS = {
        "mean_psi_per_lsv_junction": float,
        # 'stdev_psi_per_lsv_junction': float,
        "junctions_coords": str,
    }
    PSI_COORDINATE_COLUMNS = {"junctions_coords": ""}
    PSI_COLUMNS = [
        "gene_id",
        "event_type",
        "ref_exon_start",
        "ref_exon_end",
        "start",
        "end",
        "mean_psi_per_lsv_junction",
    ]
    PSI_RENAME = {"mean_psi_per_lsv_junction": "bootstrap_psi_mean"}
    return (
        df.pipe(_extract_lsv_id)
        .pipe(_df_explode, PSI_EXPLODE_COLUMNS)
        .pipe(_extract_coordinates, PSI_COORDINATE_COLUMNS)[PSI_COLUMNS]
        .rename(columns=PSI_RENAME)
    )


def extract_deltapsi(df: pd.DataFrame) -> pd.DataFrame:
    DPSI_EXPLODE_COLUMNS = {
        "mean_dpsi_per_lsv_junction": float,
        "probability_changing": float,
        "probability_non_changing": float,
        # 'stdev_psi_per_lsv_junction': float,
        "junctions_coords": str,
    }
    DPSI_COORDINATE_COLUMNS = {"junctions_coords": ""}
    DPSI_COLUMNS = [
        "gene_id",
        "event_type",
        "ref_exon_start",
        "ref_exon_end",
        "start",
        "end",
        "mean_dpsi_per_lsv_junction",
        "probability_changing",
        "probability_non_changing",
    ]
    DPSI_RENAME = {
        "mean_dpsi_per_lsv_junction": "legacy_dpsi_mean",
        "probability_changing": "legacy_probability_changing",
        "probability_non_changing": "legacy_probability_nonchanging",
    }
    return (
        df.pipe(_extract_lsv_id)
        .pipe(_df_explode, DPSI_EXPLODE_COLUMNS)
        .pipe(_extract_coordinates, DPSI_COORDINATE_COLUMNS)[DPSI_COLUMNS]
        .rename(columns=DPSI_RENAME)
    )


if __name__ == "__main__":
    import argparse
    import sys

    ANALYSIS_FUNCTIONS = {
        "psi": extract_psi,
        "deltapsi": extract_deltapsi,
    }

    parser = argparse.ArgumentParser(description="Get columns for comparison")
    parser.add_argument("legacy_tsv", type=argparse.FileType(mode="r"))
    parser.add_argument("--output-tsv", default=sys.stdout, type=argparse.FileType("w"))
    parser.add_argument(
        "--analysis-type", choices=list(ANALYSIS_FUNCTIONS.keys()), default="psi"
    )

    args = parser.parse_args()

    # load tsv
    df = pd.read_csv(args.legacy_tsv, sep="\t", comment="#")
    # later, select appropriate extraction function for analysis type
    df = ANALYSIS_FUNCTIONS[args.analysis_type](df)
    # save table
    args.output_tsv.write(f"# Produced by `{sys.argv[0]} {args.legacy_tsv.name}`\n")
    df.to_csv(args.output_tsv, sep="\t", index=False)
