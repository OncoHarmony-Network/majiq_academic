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


def assert_splicegraphs(sg1: nm.SpliceGraph, sg2: nm.SpliceGraph):
    xr.testing.assert_equal(sg1.contigs.df, sg2.contigs.df)
    xr.testing.assert_equal(sg1.genes.df, sg2.genes.df)
    xr.testing.assert_equal(sg1.exons.df, sg2.exons.df)
    xr.testing.assert_equal(sg1.introns.df, sg2.introns.df)
    xr.testing.assert_equal(sg1.junctions.df, sg2.junctions.df)


def test_gff3_command(script_runner, tmp_path):
    """Test new-majiq gff3 runs and compare to expected result"""
    path_result = str(tmp_path / "result")
    ret = script_runner.run("new-majiq", "gff3", get_path(ANNOTATED_GFF3), path_result)
    assert ret.success
    path_expected = get_path(ANNOTATED_SG)
    sg_result = nm.SpliceGraph.from_zarr(path_result)
    sg_expected = nm.SpliceGraph.from_zarr(path_expected)
    assert_splicegraphs(sg_result, sg_expected)
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
        assert_splicegraphs(sg_result, sg_expected)
    return


@pytest.mark.parametrize("min_experiments", ["1", "999"])
def test_update_vs_combine(script_runner, tmp_path, min_experiments):
    """Test majiq-build update with groups tsv vs update on sjs and combine

    We expect that running update on build groups individually and then
    combining should be equivalent to a build with them all together
    """
    grp1 = [get_sj_path(x, strandness="AUTO") for x in EXPERIMENT_GROUPS[0][1]]
    grp2 = [get_sj_path(x, strandness="AUTO") for x in EXPERIMENT_GROUPS[1][1]]
    # separate calls to majiq-build update, then majiq-build combine
    assert script_runner.run(
        "majiq-build",
        "update",
        get_path(ANNOTATED_SG),
        str(tmp_path / "grp1.sg"),
        *("--sjs", *grp1),
        "--no-simplify",
        *("--min-experiments", min_experiments),
    ).success
    assert script_runner.run(
        "majiq-build",
        "update",
        get_path(ANNOTATED_SG),
        str(tmp_path / "grp2.sg"),
        *("--sjs", *grp2),
        "--no-simplify",
        *("--min-experiments", min_experiments),
    ).success
    assert script_runner.run(
        "majiq-build",
        "combine",
        str(tmp_path / "combined.sg"),
        *("--keep-denovo", str(tmp_path / "grp1.sg"), str(tmp_path / "grp2.sg")),
    )
    # single call to majiq-build
    with open(tmp_path / "groups.tsv", mode="w") as handle:
        print("group\tsj", file=handle)
        for x in grp1:
            print(f"grp1\t{x}", file=handle)
        for x in grp2:
            print(f"grp2\t{x}", file=handle)
    assert script_runner.run(
        "majiq-build",
        "update",
        get_path(ANNOTATED_SG),
        str(tmp_path / "multigrp_build.sg"),
        *("--groups-tsv", str(tmp_path / "groups.tsv")),
        "--no-simplify",
        *("--min-experiments", min_experiments),
    ).success
    # check equivalence of splicegraphs
    assert_splicegraphs(
        nm.SpliceGraph.from_zarr(tmp_path / "combined.sg"),
        nm.SpliceGraph.from_zarr(tmp_path / "multigrp_build.sg"),
    )
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
    assert_splicegraphs(sg_result, sg_expected)
    return


@pytest.mark.parametrize("batch_group", EXPERIMENT_GROUPS)
@pytest.mark.parametrize("ignore", [None, "first", "all"])
def test_psi_coverage_command(script_runner, batch_group, ignore, tmp_path):
    """Test psi-coverage command

    Test psi-coverage command, comparing to results when not ignoring first
    experiment group events
    """
    group, names = batch_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ignore_from = list()
    if ignore:
        ignore_from = [
            "--ignore-from",
            get_path(COMBINED_SG)
            if ignore == "all"
            else get_build_path(
                EXPERIMENT_GROUPS[0][0], simplify=True, min_experiments=1
            ),
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
    if not ignore:
        path_expected = get_psicov_path(group)
        result = nm.PsiCoverage.from_zarr(path_result)
        expected = nm.PsiCoverage.from_zarr(path_expected)
        xr.testing.assert_equal(
            result.df.drop_vars(["bootstrap_psi", "bootstrap_total"]),
            expected.df.drop_vars(["bootstrap_psi", "bootstrap_total"]),
        )
    return


@pytest.mark.parametrize("batch_group", EXPERIMENT_GROUPS)
def test_sg_coverage_command(script_runner, batch_group, tmp_path):
    """Test sg-coverage command

    Smoke test of sg-coverage command
    """
    group, names = batch_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "sg-coverage",
        get_path(COMBINED_SG),
        path_result,
        *paths_sj,
    )
    assert ret.success
    return


@pytest.mark.parametrize("min_experiments", [None, 2])
def test_quantify_command(script_runner, min_experiments, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    return


# this test verifies that majiq can setup/tear down its own dask client, too
@pytest.mark.parametrize("own_cluster", [True, False])
def test_deltapsi_command(script_runner, tmp_path, dask_client, own_cluster):
    """Smoke test for deltapsi command"""
    # split up psicoverage files into two groups arbitrarily
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    half = len(paths_psicov) // 2
    dask_args = (
        ("--dask-local-directory", str(tmp_path))
        if own_cluster
        else ("--scheduler-address", dask_client.scheduler.address)
    )
    ret = script_runner.run(
        "new-majiq",
        "deltapsi",
        *("-grp1", *paths_psicov[:half]),
        *("-grp2", *paths_psicov[half:]),
        *("-n", "grp1", "grp2"),
        *("--min-experiments", "2"),
        *("--splicegraph", get_path(COMBINED_SG)),
        "--disable-progress",
        *dask_args,
    )
    assert ret.success
    return


def test_heterogen_command(script_runner, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    return


def test_psi_controls_command(script_runner, tmp_path, dask_client):
    """Smoke test for psi-controls command (all but first group)"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS[1:]]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        "new-majiq",
        "psi-controls",
        path_result,
        *paths_psicov,
        "--disable-progress",
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    return


def test_moccasin_factors_model_command(script_runner, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    xr.testing.assert_allclose(
        xr.open_zarr(path_result),
        xr.open_zarr(get_path(FACTORS_MODEL)).drop_vars(["lsv_idx", "event_size"]),
    )
    return


def test_moccasin_factors_infer_command(script_runner, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    return


def test_moccasin_coverage_model_command(script_runner, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    xr.testing.assert_allclose(
        xr.open_zarr(path_result).drop_vars(["bootstrap_model"]),
        xr.open_zarr(get_path(COVERAGE_MODEL)).drop_vars(
            ["bootstrap_model", "lsv_idx", "event_size"]
        ),
    )
    return


def test_moccasin_coverage_infer_command(script_runner, tmp_path, dask_client):
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
        *("--scheduler-address", dask_client.scheduler.address),
    )
    assert ret.success
    return
