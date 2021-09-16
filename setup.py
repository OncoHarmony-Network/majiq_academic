from setuptools import setup
from pathlib import Path

repo_path = Path(__file__).resolve().parent

setup(
    packages=[],
    install_requires=[
        f"rna_voila @ file://localhost{repo_path}/voila/",
        f"rna_majiq @ file://localhost{repo_path}/majiq/",
    ],
)
