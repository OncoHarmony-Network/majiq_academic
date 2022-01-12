"""
build_update.py

Update splicegraph with specified experiments split into build groups

Author: Joseph K Aicher
"""

import argparse
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_args import (
    BuildType,
    SJGroupsT,
    build_threshold_args,
    build_type_args,
    do_build_and_simplify,
    enable_simplify_args,
    ir_filtering_args,
    reset_simplified_args,
    simplifier_threshold_args,
)
from new_majiq.logger import get_logger

DESCRIPTION = "Update splicegraph with specified experiment groups"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "base_sg", type=ExistingResolvedPath, help="Path to base splicegraph"
    )
    parser.add_argument(
        "out_sg", type=NewResolvedPath, help="Path for output splicegraph"
    )
    experiments = parser.add_argument_group(
        "required specification of input experiments arguments"
    )
    experiments_ex = experiments.add_mutually_exclusive_group(required=True)
    experiments_ex.add_argument(
        "--groups-tsv",
        type=ExistingResolvedPath,
        metavar="TSV",
        help="Specify experiments from multiple build groups using TSV file."
        " Required columns 'group' and 'sj'."
        " One row per unique experiment."
        " `group` indicates group experiment belongs to, `sj`"
        " the path to the experiment's SJ file (from `new-majiq sj`)",
    )
    StoreSJPaths = StoreRequiredUniqueActionFactory()
    experiments_ex.add_argument(
        "--sjs",
        type=ExistingResolvedPath,
        action=StoreSJPaths,
        nargs="+",
        help="Specify experiments from a single build group directly as"
        " the unique paths to the group experiments' SJ files"
        " (use `new-majiq combine` to merge independent build groups)",
    )

    build_type_args(parser)
    ir_filtering_args(parser)
    enable_simplify_args(parser)

    reset_simplified_args(parser)
    build_threshold_args(parser)
    simplifier_threshold_args(parser, prefix="simplify-")
    resources_args(parser, use_dask=False)
    return


def get_grouped_experiments(groups_tsv: Path) -> SJGroupsT:
    """Get grouped experiments from table

    Get grouped experiments from table. Verify that there are no duplicate
    experiments (insofar as their paths) and that the paths exist.
    """
    import pandas as pd

    df = pd.read_csv(groups_tsv, sep="\t", usecols=["group", "sj"])
    if len(df.sj) == 0:
        raise ValueError(f"No input experiments specified in {groups_tsv}")
    # determine if any duplicated experiments
    duplicated_mask = df.sj.duplicated()
    if duplicated_mask.any():
        duplicated_sj = set(df.sj[duplicated_mask])
        raise ValueError(f"Requested build with repeated experiments {duplicated_sj}")
    # verify that all paths exist
    for sj_path in df.sj:
        if not Path(sj_path).exists():
            raise ValueError(f"Unable to find input experiment {sj_path}")
    return {
        group: sorted(Path(x).resolve() for x in group_sjs)
        for group, group_sjs in df.groupby("group")["sj"]
    }


def run(args: argparse.Namespace) -> None:
    simplify: bool  # will we be simplifying in the end?
    if args.simplify is None:
        simplify = args.build == BuildType.SIMPLIFY_ONLY
    else:
        if not args.simplify and args.build == BuildType.SIMPLIFY_ONLY:
            raise argparse.ArgumentError(
                None, "--no-simplify and --simplify-only are incompatible"
            )
        simplify = args.simplify

    log = get_logger()
    if args.groups_tsv:
        log.info("Loading experiment groups from %s", args.groups_tsv)
        experiments = get_grouped_experiments(args.groups_tsv)
    else:
        experiments = {"": args.sjs}
    log.info("Loading base splicegraph from %s", args.base_sg)
    sg = nm.SpliceGraph.from_zarr(args.base_sg)

    # pool for performing build, simplification
    p = ThreadPool(args.nthreads)
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
        min_experiments=args.min_experiments,
        introns_type=args.introns,
        simplify=simplify,
        reset_simplify=args.reset_simplify,
        simplify_minpsi=args.simplify_minpsi,
        simplify_minreads_annotated=args.simplify_minreads_annotated,
        simplify_minreads_denovo=args.simplify_minreads_denovo,
        simplify_minreads_ir=args.simplify_minreads_ir,
        simplify_min_experiments=args.simplify_min_experiments,
        imap_unordered_fn=p.imap_unordered,
    )
    p.close()

    log.info("Saving updated %s to %s", sg, args.out_sg)
    sg.to_zarr(args.out_sg)
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
