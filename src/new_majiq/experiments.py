"""
experiments.py

We infer experiment name information from input BAM names

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import (
    Union,
)


def bam_experiment_name(p: Union[str, Path]) -> str:
    # return file name with extension removed
    return Path(p).name.rsplit(".", 1)[0]
