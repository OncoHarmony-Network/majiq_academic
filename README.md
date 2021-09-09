# New MAJIQ!

MAJIQ is a software package for defining and quantifying RNA splicing from
RNA-seq data. This is a development branch of MAJIQ reimagining the logic of
MAJIQ v2 from the ground up using more efficient algorithms and data
structures. This was originally motivated by the need to address limitations of
how intron coverage was assessed that would be difficult to fix in the old
baseline. These changes come with the side effect of significant speed
improvements as well, along with modules that can be (or will be able to be)
used for interactive analysis (i.e. in a Jupyter notebook).

We expect this branch to eventually become MAJIQ v3.


## Installation

Requirements:

+ C++ compiler with C++17 support (e.g. gcc >= 7)
+ htslib >= 1.10
+ Python 3.8 (with a recent version of pip)
+ moccasin >= 0.26 (`new_moccasin`), originally part of majiq v3.

The installation process requires the htslib include/library directories to be
available. If they are not in `/usr/local/include` and `/usr/local/lib`, please
specify their location using the environment variables `HTSLIB_INCLUDE_DIR` and
`HTSLIB_LIBRARY_DIR`.

Install a version of moccasin since the rewrite (>= 0.26). This can currently be done
by running `pip install git+https://bitbucket.org/biociphers/moccasin@new_moccasin`.
Then, run `pip install .` in the base directory of this repo in your desired
conda/virtual environment.

If you are using conda, you can consider working with `environment.yaml`.
This is designed for development (installs in editable mode and with optional
packages) but you can edit things out further. There will be a separate
environment in the future that can be used with Snakemake to automatically
create an appropriate environment if SSH keys and environment variables are set
up appropriately.

### For development

Install `pre-commit` using pip (or use `environment.yaml`), then run
`pre-commit install`.
New changes will automatically be formatted/checked for potential issues.
It is occasionally worth checking all files with these checks -- this can be
done by running `pre-commit run --all-files`.


## Running

TODO
