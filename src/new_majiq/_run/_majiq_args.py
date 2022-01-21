"""
majiq_args.py

Helper classes/functions for checking valid arguments directly at command-line
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Callable, Optional, Set, Union

import new_majiq as nm


def ExistingResolvedPath(x: Any) -> Path:
    try:
        p = Path(x).resolve()
    except TypeError:
        raise argparse.ArgumentTypeError(f"Input {x} is not pathlike")
    if not p.exists():
        raise argparse.ArgumentTypeError(f"{p} does not exist")
    return p


def NewResolvedPath(x: Any) -> Path:
    try:
        p = Path(x).resolve()
    except TypeError:
        raise argparse.ArgumentTypeError(f"Input {x} is not pathlike")
    if p.exists():
        raise argparse.ArgumentTypeError(f"{p} already exists")
    return p


def check_characters_factory(
    valid_characters: Set[str],
    valid_description: Optional[str] = None,
) -> Callable[[Any], str]:
    """Create function to raise argparse type error or return string

    Create function to raise argparse.ArgumentTypeError if input contains
    invalid characters, otherwise, return str(x)

    Parameters
    ----------
    valid_characters: Set[str]
        characters that are acceptable for this restricted string
    valid_description: Optional[str]
        description of acceptable characters in input
    """

    def check_characters(x: Any) -> str:
        x = str(x)
        set_x = set(x)
        if not set_x.issubset(valid_characters):
            invalid_chars = set_x.difference(valid_characters)
            prefix = (
                f"may only contain {valid_description}: " if valid_description else ""
            )
            raise argparse.ArgumentTypeError(
                f"{prefix}{x} contains invalid characters {invalid_chars}"
            )
        return x

    return check_characters


def check_range_factory(
    cast_fn: Callable[[Any], Union[int, float]],
    minval: Union[int, float] = float("-inf"),
    maxval: Union[int, float] = float("inf"),
    mininclude: bool = False,
    maxinclude: bool = False,
):
    """Return argparse type function to casted type in specified range"""
    not_too_small = (lambda x: x >= minval) if mininclude else (lambda x: x > minval)
    not_too_big = (lambda x: x <= maxval) if maxinclude else (lambda x: x < maxval)
    invalid_message = (
        "parameter must be in range"
        f" {'[' if mininclude else '('}{minval}, {maxval}{']' if maxinclude else ')'}"
    )

    def return_function(value):
        ivalue = cast_fn(value)
        if not_too_small(ivalue) and not_too_big(ivalue):
            return ivalue
        else:
            raise argparse.ArgumentTypeError(
                f"{value} outside of valid range; {invalid_message}"
            )

    return return_function


def check_nonnegative_factory(
    cast_fn: Callable[[Any], Union[int, float]], reject_zero: bool
):
    """Returns argparse type function to casted type, rejecting negative values"""
    return check_range_factory(cast_fn, minval=0, mininclude=not reject_zero)


def StoreRequiredUniqueActionFactory():
    """Return class that pools shared values (do not allow overlap)"""

    class StoreRequiredUniqueAction(argparse.Action):

        # class variable with shared values known to all instances
        shared_values = {}  # Dict[str, Set[str]]

        @staticmethod
        def repeated_values(values):
            """Returns repeated values (that occur more than twice)

            Returns
            -------
            Set[str]
            """
            unique_values = set(values)  # set of unique values
            repeated_values = set()  # set of values that were repeated
            if isinstance(values, list) and len(values) != len(unique_values):
                # values are not all unique...
                repeated_values = set()  # determine which values were repeated
                # remove unique values (if done twice, added to repeated_values)
                for x in values:
                    try:
                        unique_values.remove(x)
                    except KeyError:
                        repeated_values.add(x)
            return repeated_values

        @classmethod
        def overlapping_values(cls, dest, values):
            """Updates shared_values and returns overlapping values with other
            parameters

            Parameters
            ----------
            dest: str
                Key to ignore in shared_values
            values:
                Values to compare with other shared values

            Returns
            -------
            Dict[str, Set[str]]
                Keys for other parameter names, values correspond to overlaps
                between them
            """
            unique_values = set(values)
            overlapping_values = {
                other_key: other_values & unique_values
                for other_key, other_values in cls.shared_values.items()
                if other_key != dest and other_values & unique_values
            }
            # update shared values
            cls.shared_values[dest] = unique_values
            # remove any keys with zero overlaps
            return {k: v for k, v in overlapping_values.items() if v}

        def __call__(self, parser, namespace, values, option_string=None):
            """Check that values are unique if list, not in shared_values"""
            # no repeated or overlapping values
            repeated_values = self.repeated_values(values)
            if repeated_values:
                raise argparse.ArgumentError(
                    self, f"values must be unique (repeated values: {repeated_values})"
                )
            overlapping_values = self.overlapping_values(self.dest, values)
            if overlapping_values:
                raise argparse.ArgumentError(
                    self,
                    "values must not overlap shared parameters"
                    f" (overlaps: {overlapping_values})",
                )
            # no non-unique or overlapping values, so save updated values
            setattr(namespace, self.dest, values)

    return StoreRequiredUniqueAction


def _quantify_shared_args(parser: argparse.ArgumentParser) -> None:
    """shared arguments for quantify_?(no)comparison_args"""
    parser.add_argument(
        "--splicegraph",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="If specified, annotate quantifications with splicegraph information",
    )
    parser.add_argument(
        "--output-tsv",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )


def quantify_nocomparison_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to PsiCoverage files to quantify",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        default=None,
        help="If specified, treat samples as replicates and quantify combined"
        " coverage for events that passed in at least min experiments"
        " (proportion of experiments if < 1)."
        " Default is to quantify each sample independently.",
    )
    _quantify_shared_args(parser)
    return


def quantify_comparison_args(parser: argparse.ArgumentParser) -> None:
    comparison_req = parser.add_argument_group(
        "required quantification group arguments"
    )
    StorePSICovPaths = StoreRequiredUniqueActionFactory()
    comparison_req.add_argument(
        "-psi1",
        "-grp1",
        "-g1",
        dest="psi1",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        required=True,
        help="Paths to PsiCoverage files for experiments contributing to psi1",
    )
    comparison_req.add_argument(
        "-psi2",
        "-grp2",
        "-g2",
        dest="psi2",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        required=True,
        help="Paths to PsiCoverage files for experiments contributing to psi2",
    )
    StoreGroupNames = StoreRequiredUniqueActionFactory()
    check_group_chars = check_characters_factory(
        nm.constants.ALLOWED_GROUP_NAME_CHARS, "alphanumeric or underscore characters"
    )
    comparison_req.add_argument(
        "-n",
        "--names",
        dest="names",
        nargs=2,
        metavar=("NAME1", "NAME2"),
        required=True,
        action=StoreGroupNames,
        type=check_group_chars,
        help="The names that identify the groups being compared.",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. This specifies the fraction"
        " (value < 1) or absolute number (value >= 1) of experiments that must"
        " pass individually in each group for an LSV to be quantified"
        " (default: %(default)s)",
    )
    _quantify_shared_args(parser)
    parser.add_argument(
        "--output-voila",
        metavar="ZARR",
        type=NewResolvedPath,
        default=None,
        help="Path for output voila view file (default: none)",
    )
    chunks_args(parser, nm.constants.DEFAULT_COVERAGE_CHUNKS)
    return


def chunks_args(parser: argparse.ArgumentParser, chunksize: int) -> None:
    chunks = parser.add_argument_group("output chunks arguments")
    chunks.add_argument(
        "--chunksize",
        metavar="CHUNKS",
        type=check_nonnegative_factory(int, True),
        default=chunksize,
        help="Chunk output per prefix and per this many introns/junctions at"
        " a time (default: %(default)s)",
    )


def resources_args(parser: argparse.ArgumentParser, use_dask: bool = False) -> None:
    resources = parser.add_argument_group("threads/resources arguments")
    if use_dask:
        resources.add_argument(
            "--nthreads",
            "-j",
            metavar="N",
            type=check_nonnegative_factory(int, True),
            default=nm.constants.DEFAULT_QUANTIFY_NTHREADS,
            help="Number of threads to perform work in chunks for Dask scheduler"
            " (default: %(default)s)",
        )
        resources.add_argument(
            "--memory-limit",
            metavar="MEM",
            type=str,
            default="auto",
            help="Memory limit to pass to dask cluster (default: %(default)s)",
        )
        resources.add_argument(
            "--dask-local-directory",
            metavar="path",
            type=Path,
            default=None,
            help="Path for dask-worker-space (default: current directory)",
        )
        progress_bar = resources.add_mutually_exclusive_group()
        progress_bar.add_argument(
            "--show-progress",
            dest="show_progress",
            action="store_true",
            default=True,
            help="Show progress bar for Dask computations"
            " (default: show_progress=%(default)s)",
        )
        progress_bar.add_argument(
            "--disable-progress",
            dest="show_progress",
            action="store_false",
            default=True,
            help="Disable progress bar for Dask computations"
            " (default: show_progress=%(default)s)",
        )
        resources.add_argument(
            "--scheduler-address",
            metavar="url",
            type=str,
            default=None,
            help="Use resources from existing Dask cluster at specified URL."
            " This ignores nthreads, memory_limit, and dask_local_directory."
            " Default is to create a new local cluster with those resources"
            " for the duration of the script."
            " If used, expects cluster to remain running and accessible while"
            " the script is running.",
        )
        parser.set_defaults(use_dask=True)
    else:
        resources.add_argument(
            "--nthreads",
            "-j",
            metavar="N",
            type=check_nonnegative_factory(int, True),
            default=nm.constants.DEFAULT_BAM_NTHREADS,
            help="Number of threads used for simultaneous processing of multiple"
            " input files (default: %(default)s)",
        )
