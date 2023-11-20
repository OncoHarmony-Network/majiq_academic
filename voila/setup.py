import os
from setuptools import setup


def files_recursive(directory):
    return [
        os.path.join("..", path, filename)
        for path, directories, filenames in os.walk(directory)
        for filename in filenames
    ]


setup(
    package_data={
        "rna_voila": [
            "api/model.sql",
            *files_recursive("rna_voila/view/templates"),
            *files_recursive("rna_voila/view/static"),
        ],
    },
    use_scm_version={
        "root": "..",
        "relative_to": __file__,
        "fallback_version": "2.5.0+scmfallback",
    },
)
