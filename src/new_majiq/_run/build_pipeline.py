"""
build_pipeline.py

Run commands of majiq-build in sequence:
+ Load annotated splicegraph (majiq-build gff3)
+ Load coverage from BAM files (majiq-build sj)
+ Update annotated splicegraph (majiq-build update)
+ Obtain coverage over LSVs (majiq-build psi-coverage)

In the future, we might update coverage to also have splicegraph coverage in
the same file
"""

import argparse
import multiprocessing.dummy as mp
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import (
    SJGroupsT,
    build_threshold_args,
    build_type_args,
    do_build_and_simplify,
    enable_simplify_args,
    gff3_parsing_args,
    gff3_parsing_pipeline,
    ir_filtering_args,
    lsv_coverage_args,
    quantifiability_threshold_args,
    simplifier_threshold_args,
    sj_strandness_args,
)
from new_majiq.experiments import bam_experiment_name
from new_majiq.logger import get_logger

DESCRIPTION = (
    "majiq-build pipeline to build splicegraph and coverage from annotations"
    " and RNA-seq experiments"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "gff3",
        type=ExistingResolvedPath,
        help="Path to GFF3 file (uncompressed or gzipped) with annotated"
        " transcripts for initial splicegraph",
    )
    parser.add_argument(
        "experiments_tsv",
        type=ExistingResolvedPath,
        help="Path to TSV file indicating experiments from one or more build"
        " groups."
        " Required column 'path', indicating path to input file for each unique"
        " experiment."
        " Column `group` indicates group experiment belongs to; `path`"
        " indicates the path to the input file used for the experiment."
        " Optional columns 'group' and 'strandness'."
        " Column `group` indicates group experiment belongs to, by default all"
        " experiments are in a single build group."
        " Column `strandness` overrides value from --strandness."
        " It is assumed that previously loaded SJExperiment have extension"
        " `sj`, `zarr`, or `zip` (and alignments do not).",
    )
    parser.add_argument(
        "output_dir",
        type=NewResolvedPath,
        help="New path for output directory with files produced by the pipeline",
    )
    build_type_args(parser, allow_simplify_only=False)
    ir_filtering_args(parser)
    enable_simplify_args(parser)
    build_threshold_args(parser)
    simplifier_threshold_args(parser, prefix="simplify-")
    quantifiability_threshold_args(parser, prefix="quantify-")
    lsv_coverage_args(parser)
    sj_strandness_args(parser)
    gff3_parsing_args(parser, prefix="gff3-")
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    import pandas as pd

    # load args.experiments_tsv
    df_experiments = pd.read_csv(
        args.experiments_tsv,
        sep="\t",
        dtype={
            "path": str,
            "group": str,
            "strandness": str,
            "is_sj": bool,
        },
    )
    if "path" not in df_experiments.columns:
        raise ValueError("Experiments table does not have 'path' column")
    elif len(set(df_experiments["path"])) < len(df_experiments):
        nonunique_paths = (
            df_experiments["path"]
            .value_counts()
            .pipe(lambda x: x[x > 1])
            .index.tolist()
        )
        raise ValueError(f"Input paths are not unique ({nonunique_paths = })")
    else:
        df_experiments["path"] = df_experiments["path"].apply(ExistingResolvedPath)
    if "group" not in df_experiments.columns:
        df_experiments["group"] = "experiments"  # all are same group, trivial name
    if "strandness" not in df_experiments.columns:
        df_experiments["strandness"] = args.strandness
    else:
        df_experiments.fillna({"strandness": args.strandness}, inplace=True)
        # all should be uppercase
        df_experiments["strandness"] = df_experiments["strandness"].str.upper()
        valid_strandness = {"AUTO", *nm.ExperimentStrandness.__entries.keys()}
        if invalid_strandness := set(df_experiments["strandness"]) - valid_strandness:
            raise ValueError(
                "Invalid values for strandness in experiments_tsv"
                f" ({invalid_strandness = }, {valid_strandness = })"
            )
    df_experiments["is_sj"] = df_experiments["path"].apply(
        lambda x: x.suffix.upper() in {".SJ", ".ZIP", ".ZARR"}
    )
    df_experiments["sj_path"] = [
        path if is_sj else args.output_dir / f"{bam_experiment_name(path)}.sj"
        for path, is_sj in zip(df_experiments["path"], df_experiments["is_sj"])
    ]

    log = get_logger()
    log.info(f"Creating directory for output files at {args.output_dir}")
    args.output_dir.mkdir()

    # load base splicegraph
    log.info(
        "Load annotated splicegraph from GFF3"
        f" (`majiq-build gff3 {str(args.gff3)} <memory>`)"
    )
    sg: nm.SpliceGraph = gff3_parsing_pipeline(
        args.gff3,
        features_default=args.gff3_features_default,
        features_tsv=args.gff3_features_tsv,
        types_genes=args.gff3_types_genes,
        types_transcripts=args.gff3_types_transcripts,
        types_exons=args.gff3_types_exons,
        types_silent=args.gff3_types_silent,
        types_hard_skip=args.gff3_types_hard_skip,
    )

    # get SJ files that we don't already have
    log.info("Preparing SJ coverage for each experiment")
    sj: Optional[nm.SJExperiment] = None
    df_experiments_bam = df_experiments.loc[~df_experiments["is_sj"]]
    for idx, (_, experiment_info) in enumerate(df_experiments_bam.iterrows(), 1):
        log.info(
            f"(`majiq-build sj {experiment_info['path']} <memory>"
            f" {experiment_info['sj_path']}`)"
        )
        sj = nm.SJExperiment.from_bam(
            experiment_info["path"],
            sg,
            experiment_info["strandness"],
            nthreads=args.nthreads,
            auto_minreads=args.auto_minreads,
            auto_minjunctions=args.auto_minjunctions,
            auto_mediantolerance=args.auto_mediantolerance,
        )
        log.info(
            f"Saving SJExperiment from {experiment_info['path']}"
            f" to {experiment_info['sj_path']} ({idx} / {len(df_experiments_bam)})"
        )
        sj.to_zarr(experiment_info["sj_path"])
    sj = None  # unload SJ coverage from memory

    # create pool for multithreading
    p = mp.Pool(args.nthreads)

    output_splicegraph = args.output_dir / "splicegraph.zarr"
    log.info(
        f"(`majiq-build update <memory> {output_splicegraph} --groups-tsv <memory>`)"
    )
    experiments: SJGroupsT = {
        group: sj_paths.tolist()
        for group, sj_paths in df_experiments.groupby("group")["sj_path"]
    }
    sg = do_build_and_simplify(
        sg,
        experiments,
        build=args.build,
        minreads=args.minreads,
        mindenovo=args.mindenovo,
        minpos=args.minpos,
        max_pctbins=args.max_pctbins,
        junction_acceptance_probability=args.junction_acceptance_probability,
        intron_acceptance_probability=args.intron_acceptance_probability,
        simplify=bool(args.simplify),
        simplify_minpsi=args.simplify_minpsi,
        simplify_minreads_annotated=args.simplify_minreads_annotated,
        simplify_minreads_denovo=args.simplify_minreads_denovo,
        simplify_minreads_ir=args.simplify_minreads_ir,
        simplify_min_experiments=args.simplify_min_experiments,
        imap_unordered_fn=p.imap_unordered,
    )
    log.info(f"Saving updated splicegraph to {output_splicegraph}")
    sg.to_zarr(output_splicegraph)

    log.info(f"Defining LSVs for coverage ({args.select_lsvs})")
    lsvs = sg.exon_connections.lsvs(args.select_lsvs)
    if args.ignore_from is not None:
        log.info(f"Ignoring LSVs also found in {args.ignore_from}")
        lsvs = lsvs[
            lsvs.unique_events_mask(
                nm.SpliceGraph.from_zarr(
                    args.ignore_from, genes=sg.genes
                ).exon_connections.lsvs(args.select_lsvs)
            ).unique_events_mask
        ]

    nm.rng_resize(args.nthreads)
    for group_idx, (group, sj_paths) in enumerate(experiments.items(), 1):
        output_psicov = args.output_dir / f"{group}.psicov"
        sj_paths_str: str = str(sj_paths[0])
        if len(sj_paths) > 2:
            sj_paths_str += " {...}"
        if len(sj_paths) > 1:
            sj_paths_str += f" {str(sj_paths[-1])}"
        log.info(f"Saving PsiCoverage for {group} ({group_idx} / {len(experiments)})")
        log.info(
            f" (`majiq-build psi-coverage {output_splicegraph} {output_psicov}"
            f" {sj_paths_str}`)"
        )
        nm.PsiCoverage.convert_sj_batch(
            sj_paths,
            lsvs,
            output_psicov,
            minreads=args.quantify_minreads,
            minbins=args.quantify_minbins,
            num_bootstraps=args.num_bootstraps,
            pvalue_threshold=args.stack_pvalue_threshold,
            imap_unordered_fn=p.imap_unordered,
        )
    return


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args(sys_args)
    run(args)
    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


if __name__ == "__main__":
    main()
