"""
conftest.py

This is automatically used by all tests in this directory for fixture definition

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Tuple

import pandas as pd
import pytest
from dask.distributed import Client


@pytest.fixture(scope="session")
def dask_client(tmpdir_factory):
    # not multiprocess, session-duration local directory for workers
    worker_space = tmpdir_factory.mktemp("worker-space")
    client = Client(n_workers=1, local_directory=worker_space)
    # yield client to different tests
    yield client
    # clean up after yield
    client.shutdown()
    del client
    return


@pytest.fixture(scope="session")
def het_vignette_data(tmpdir_factory, request) -> Tuple[Path, pd.DataFrame]:
    vignette_dir = Path(__file__).parent / "data/het_vignette"
    if not vignette_dir.exists():
        import shutil
        import urllib.request
        import zipfile

        capmanager = request.config.pluginmanager.getplugin("capturemanager")
        with capmanager.global_and_fixture_disabled():
            tmp_path = Path(tmpdir_factory.mktemp("archive_download"))
            print("Downloading MAJIQ HET vignette data for testing")
            zip_path = Path(
                urllib.request.urlretrieve(
                    "http://majiq.biociphers.org/data/majiq_het_vignette.zip",
                    tmp_path / "het_vignette.zip",
                )[0]
            )
            print("Extracting vignette data")
            extract_dir = tmp_path / "extract_dir"
            with zipfile.ZipFile(zip_path, mode="r") as f:
                f.extractall(extract_dir)
            zip_path.unlink()  # delete zip archive
            print("Deleting unnecessary data from vignette archive")
            (extract_dir / "settings.ini").unlink()
            (extract_dir / "vignette.ipynb").unlink()
            shutil.rmtree(extract_dir / "images")
            shutil.rmtree(extract_dir / "output")
            print(f"Moving data for testing to {vignette_dir}")
            extract_dir.rename(vignette_dir)
            print("Done preparing vignette data for testing")
    # prepare information that tests need
    gff3_path = vignette_dir / "gff3/ensembl.Homo_sapiens.GRCh38.94.subset.gff3.gz"
    df_experiments = (
        pd.read_csv(vignette_dir / "build-info.tsv", sep="\t")
        .assign(
            comparison_group=lambda df: df.apply(
                lambda x: f"{x.cell_line}-{x.knockdown}", axis=1
            ),
            bam_path=lambda df: df.prefix.apply(
                lambda x: str(vignette_dir / f"bam/{x}.subset.bam")
            ),
            prefix=lambda df: df.bam_path.apply(
                lambda x: Path(x).name.rsplit(".", 1)[0]
            ),
        )
        .rename(columns={"encode_accession": "replicate_group"})[
            ["comparison_group", "replicate_group", "prefix", "bam_path"]
        ]
    )
    return gff3_path, df_experiments
