"""
test_scripts.py

Test of console scripts (mostly smoke tests, but also checking against
saved/expected output in select cases)

Depends on pytest_console_scripts for script_runner fixture

Author: Joseph K Aicher
"""

from pathlib import Path

import pytest
import xarray as xr

import new_majiq as nm
from new_majiq.experiments import bam_experiment_name

ANNOTATED_GFF3 = "data/gff3/SRSF4.gff3"
ANNOTATED_SG = "data/sg/annotated.sg.zip"
EXPERIMENT_NAMES = sorted(
    bam_experiment_name(x) for x in Path(__file__).parent.glob("data/bam/*.bam")
)
EXPERIMENT_GROUPS = [
    ("DUM", [x for x in EXPERIMENT_NAMES if "DUM" in x]),
    ("TIZ", [x for x in EXPERIMENT_NAMES if "TIZ" in x]),
    ("JBU", [x for x in EXPERIMENT_NAMES if "JBU" in x]),
    ("AOV", [x for x in EXPERIMENT_NAMES if "AOV" in x]),
]
COMBINED_SG = "data/sg/combined.sg.zip"
FACTORS_MODEL = "data/factors/factor_model.zarr.zip"
FACTORS_TSV = "data/factors/factors.tsv"
FACTORS_TSV_CONFOUNDERS = ["confounding"]
COVERAGE_MODEL = "data/factors/coverage_model.zarr.zip"


def get_path(p) -> str:
    """Path to p relative to test script"""
    return str(Path(__file__).parent / p)


def get_bam_path(name):
    """Path to BAM files"""
    return get_path(f"data/bam/{name}.bam")


def get_sj_path(name, strandness):
    """Path to SJ files as input/expected from BAM"""
    return get_path(f"data/sj/{name}.{strandness}.sj.zip")


def get_build_path(group, simplify, min_experiments):
    """Path to splicegraph files as input/expected from build command"""
    return get_path(f"data/build/{group}.{simplify}.{min_experiments}.sg.zip")


def get_psicov_path(group):
    """Path to psicoverage files as input/expected from psi-coverage command"""
    return get_path(f"data/psicov/{group}.psicov.zip")


def test_gff3_command(script_runner, tmp_path):
    """Test new-majiq gff3 runs and compare to expected result"""
    path_result = str(tmp_path / "result")
    ret = script_runner.run("new-majiq", "gff3", get_path(ANNOTATED_GFF3), path_result)
    assert ret.success
    path_expected = get_path(ANNOTATED_SG)
    sg_result = nm.SpliceGraph.from_zarr(path_result)
    sg_expected = nm.SpliceGraph.from_zarr(path_expected)
    xr.testing.assert_equal(sg_result.contigs.df, sg_expected.contigs.df)
    xr.testing.assert_equal(sg_result.genes.df, sg_expected.genes.df)
    xr.testing.assert_equal(sg_result.exons.df, sg_expected.exons.df)
    xr.testing.assert_equal(sg_result.introns.df, sg_expected.introns.df)
    xr.testing.assert_equal(sg_result.junctions.df, sg_expected.junctions.df)
    return


@pytest.mark.parametrize("name", EXPERIMENT_NAMES)
@pytest.mark.parametrize("strandness", ["AUTO", "NONE", "FORWARD", "REVERSE"])
def test_sj_command(script_runner, name, strandness, tmp_path):
    """Test new-majiq sj command. Check results if found in data (AUTO)"""
    path_bam = get_bam_path(name)
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq", "sj", path_bam, get_path(ANNOTATED_SG), path_result
    )
    assert ret.success
    path_expected = get_sj_path(name, strandness=strandness)
    if Path(path_expected).exists():
        sj_result = nm.SJExperiment.from_zarr(path_result)
        sj_expected = nm.SJExperiment.from_zarr(path_expected)
        xr.testing.assert_equal(sj_result.junctions.df, sj_expected.junctions.df)
        xr.testing.assert_equal(sj_result.introns.df, sj_expected.introns.df)
    return


@pytest.mark.parametrize("build_group", EXPERIMENT_GROUPS)
@pytest.mark.parametrize("simplify", [False, True])
@pytest.mark.parametrize("min_experiments", [0.5, 1])
def test_build_command(script_runner, build_group, simplify, min_experiments, tmp_path):
    """Test new-majiq build command

    Test new-majiq build command. Check results if found in data
    (when simplify=True, min-experiments=1)
    """
    group, names = build_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "build",
        get_path(ANNOTATED_SG),
        path_result,
        "--sjs",
        *paths_sj,
        "--simplify" if simplify else "--no-simplify",
        "--min-experiments",
        str(min_experiments),
    )
    assert ret.success
    path_expected = get_build_path(
        group, simplify=simplify, min_experiments=min_experiments
    )
    if Path(path_expected).exists():
        sg_result = nm.SpliceGraph.from_zarr(path_result)
        sg_expected = nm.SpliceGraph.from_zarr(path_expected)
        xr.testing.assert_equal(sg_result.contigs.df, sg_expected.contigs.df)
        xr.testing.assert_equal(sg_result.genes.df, sg_expected.genes.df)
        xr.testing.assert_equal(sg_result.exons.df, sg_expected.exons.df)
        xr.testing.assert_equal(sg_result.introns.df, sg_expected.introns.df)
        xr.testing.assert_equal(sg_result.junctions.df, sg_expected.junctions.df)
    return


def test_combine_command(script_runner, tmp_path):
    """Test new-majiq combine command

    Test new-majiq combine command. Set first experiment group as annotated.
    Compare result to expected.
    """
    make_annotated = [
        get_build_path(group, simplify=True, min_experiments=1)
        for group, _ in EXPERIMENT_GROUPS[:1]
    ]
    keep_denovo = [
        get_build_path(group, simplify=True, min_experiments=1)
        for group, _ in EXPERIMENT_GROUPS[1:]
    ]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "combine",
        path_result,
        "--make-annotated",
        *make_annotated,
        "--keep-denovo",
        *keep_denovo,
    )
    assert ret.success
    path_expected = get_path(COMBINED_SG)
    sg_result = nm.SpliceGraph.from_zarr(path_result)
    sg_expected = nm.SpliceGraph.from_zarr(path_expected)
    xr.testing.assert_equal(sg_result.contigs.df, sg_expected.contigs.df)
    xr.testing.assert_equal(sg_result.genes.df, sg_expected.genes.df)
    xr.testing.assert_equal(sg_result.exons.df, sg_expected.exons.df)
    xr.testing.assert_equal(sg_result.introns.df, sg_expected.introns.df)
    xr.testing.assert_equal(sg_result.junctions.df, sg_expected.junctions.df)
    return


@pytest.mark.parametrize("batch_group", EXPERIMENT_GROUPS)
@pytest.mark.parametrize("ignore_first", [False, True])
def test_psi_coverage_command(script_runner, batch_group, ignore_first, tmp_path):
    """Test psi-coverage command

    Test psi-coverage command, comparing to results when not ignoring first
    experiment group events
    """
    group, names = batch_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ignore_from = list()
    if ignore_first:
        ignore_from = [
            "--ignore-from",
            get_build_path(EXPERIMENT_GROUPS[0][0], simplify=True, min_experiments=1),
        ]
    ret = script_runner.run(
        "new-majiq",
        "psi-coverage",
        get_path(COMBINED_SG),
        path_result,
        *paths_sj,
        *ignore_from,
    )
    assert ret.success
    if not ignore_first:
        path_expected = get_psicov_path(group)
        result = nm.PsiCoverage.from_zarr(path_result)
        expected = nm.PsiCoverage.from_zarr(path_expected)
        xr.testing.assert_equal(
            result.df.drop_vars(["bootstrap_psi", "bootstrap_total"]),
            expected.df.drop_vars(["bootstrap_psi", "bootstrap_total"]),
        )
    return


@pytest.mark.parametrize("min_experiments", [None, 2])
def test_quantify_command(script_runner, min_experiments, tmp_path):
    """Smoke test for quantify command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    flag = []
    if min_experiments:
        flag = ["--min-experiments", str(min_experiments)]
    ret = script_runner.run(
        "new-majiq",
        "quantify",
        *paths_psicov,
        *flag,
        *("--splicegraph", get_path(COMBINED_SG)),
        *("--quantiles", "0.1", "0.9"),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return


def test_deltapsi_command(script_runner, tmp_path):
    """Smoke test for deltapsi command"""
    # split up psicoverage files into two groups arbitrarily
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    half = len(paths_psicov) // 2
    ret = script_runner.run(
        "new-majiq",
        "deltapsi",
        *("-grp1", *paths_psicov[:half]),
        *("-grp2", *paths_psicov[half:]),
        *("-n", "grp1", "grp2"),
        *("--min-experiments", "2"),
        *("--splicegraph", get_path(COMBINED_SG)),
        *("--use-posterior", "both"),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return


def test_heterogen_command(script_runner, tmp_path):
    """Smoke test for heterogen command"""
    # split up psicoverage files into two groups arbitrarily
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    half = len(paths_psicov) // 2
    ret = script_runner.run(
        "new-majiq",
        "heterogen",
        *("-grp1", *paths_psicov[:half]),
        *("-grp2", *paths_psicov[half:]),
        *("-n", "grp1", "grp2"),
        *("--min-experiments", "2"),
        *("--splicegraph", get_path(COMBINED_SG)),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return


def test_psi_controls_command(script_runner, tmp_path):
    """Smoke test for psi-controls command (all but first group)"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS[1:]]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "psi-controls",
        path_result,
        *paths_psicov,
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return


def test_moccasin_factors_model_command(script_runner, tmp_path):
    """Smoke test for moccasin-factors-model command, check close to expected"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "moccasin-factors-model",
        path_result,
        *paths_psicov,
        "--intercept-only",
        *("--ruv-max-new-factors", "1"),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    xr.testing.assert_allclose(
        xr.open_zarr(path_result), xr.open_zarr(get_path(FACTORS_MODEL))
    )
    return


def test_moccasin_factors_infer_command(script_runner, tmp_path):
    """Smoke test for moccasin-factors-infer command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "moccasin-factors-infer",
        get_path(FACTORS_MODEL),
        path_result,
        *paths_psicov,
        "--intercept-only",
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return


def test_moccasin_coverage_model_command(script_runner, tmp_path):
    """Smoke test for moccasin-coverage-model command, check close to expected"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "moccasin-coverage-model",
        path_result,
        *paths_psicov,
        *("--factors-tsv", get_path(FACTORS_TSV)),
        *("--confounding", *FACTORS_TSV_CONFOUNDERS),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    xr.testing.assert_allclose(
        xr.open_zarr(path_result).drop_vars(["bootstrap_model"]),
        xr.open_zarr(get_path(COVERAGE_MODEL)).drop_vars(["bootstrap_model"]),
    )
    return


def test_moccasin_coverage_infer_command(script_runner, tmp_path):
    """Smoke test for moccasin-coverage-infer command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "moccasin-coverage-infer",
        get_path(COVERAGE_MODEL),
        path_result,
        *paths_psicov,
        *("--factors-tsv", get_path(FACTORS_TSV)),
        *("--confounding", *FACTORS_TSV_CONFOUNDERS),
        "--disable-progress",
        *("--dask-local-directory", str(tmp_path)),
    )
    assert ret.success
    return