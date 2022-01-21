"""
conftest.py

This is automatically used by all tests in this directory for fixture definition

Author: Joseph K Aicher
"""

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
