# New MAJIQ!

MAJIQ uses the results of RNA-seq experiments to analyze alternative splicing.
This is a development branch of MAJIQ reimagining the logic of MAJIQ v2 from
the ground up using more efficient algorithms and data structures.
This was originally motivated by the need to address fundamental limitations of
how intron coverage was measured in MAJIQ v2.
Additional influence since has been by the MAJIQ-CLIN project.
These changes come along with significant speed improvements along with modules
that can be used for more interactive analysis (i.e. in a Jupyter notebook).

We expect this branch to eventually become MAJIQ v3.


# Installation

## Requirements

MAJIQ requires:

+ [HTSLib][htslib-src] (1.10 or later).
+ C++ compiler with C++17 (e.g. gcc >= 7)
+ Python (3.8 or later).
+ setuptools (45 or later).
+ moccasin >= 0.26 (`new_moccasin`), originally part of majiq v3.


## HTSLib installation

MAJIQ requires [HTSLib][htslib-src] (>= 1.10) to be installed.
This may already available in a shared directory of your environment,
otherwise, you can install HTSLib from source.
For example, to install htslib-1.13 to `~/install/htslib-1.13`, you could run:

```bash
# download htslib 1.13 archive to current directory
curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xf htslib-1.13.tar.bz2  # extract archive
cd htslib-1.13  # change into htslib source directory
# configure, make, and install htslib to ~/install/htslib-1.13
./configure --prefix=$HOME/install/htslib-1.13
make
make install
```

See [`INSTALL` from HTSLib source][htslib-install] for detailed installation
instructions for HTSLib.

   [htslib-src]: https://github.com/samtools/htslib/releases
   [htslib-install]: https://raw.githubusercontent.com/samtools/htslib/1.13/INSTALL


## Pip installation

Install MAJIQ using pip.
MAJIQ setup needs to know where HTSLib is installed.
By default, it assumes that the library and include directories are in the Unix
default locations: `/usr/local/lib`, `/usr/local/include`.
If this is not the case, please set the environment variables
`HTSLIB_LIBRARY_DIR` and `HTSLIB_INCLUDE_DIR` to the actual locations.
Install a version of moccasin since the rewrite (>= 0.26). This can currently
be done by running
`pip install git+https://bitbucket.org/biociphers/moccasin@new_moccasin`.
Then install this package from base directory of this repository:
`pip install .`.

For example, with previous instructions, from same directory as README:

```bash
# change to where library/include directories are installed
export HTSLIB_LIBRARY_DIR=$HOME/install/htslib-1.13/lib
export HTSLIB_INCLUDE_DIR=$HOME/install/htslib-1.13/include
# install moccasin
pip install git+https://bitbucket.org/biociphers/moccasin@new_moccasin
# NOTE: from same directory as this README
pip install .  # install both majiq and voila
```


## Conda installation

If you are using conda, you can consider working with `environment.yaml`
(`conda env create -n {DESIRED_ENV_NAME} -f environment.yaml`).
This is designed for development (installs in editable mode and with optional
packages) but you can edit things out further. There will be a separate
environment in the future that can be used with Snakemake to automatically
create an appropriate environment if SSH keys and environment variables are set
up appropriately.


## For development

Install `pre-commit` using pip (or use `environment.yaml`), then run
`pre-commit install`.
New changes will automatically be formatted/checked for potential issues.
It is occasionally worth checking all files with these checks -- this can be
done by running `pre-commit run --all-files`.


# Usage

## MAJIQ Builder

New MAJIQ/v3 is organized differently than MAJIQ v2.
MAJIQ v2 was organized into MAJIQ builder vs MAJIQ quantifiers.
The MAJIQ v2 builder performed several tasks in sequence, which are separated in
`new-majiq`:

1. Load/infer annotated splicegraph from GFF3 annotations (`new-majiq gff3`).
2. Parse experiment junction/intron coverage from aligned BAMs (`new-majiq sj`).
3. Use coverage at junctions and introns to infer updated splicegraph
   (`new-majiq build`/`new-majiq combine`).
4. Get coverage per LSV from splicegraph/SJ files (`new-majiq psi-coverage`).


For example, if we wanted to analyze splicing in `bam/example{1,2}.bam`
relative to annotations in `gff3/example.gff3`, we would run:

```bash
# load example.gff3
# create annotated splicegraph file at results/annotated.splicegraph
new-majiq gff3 gff3/example.gff3 results/annotated.splicegraph

# load example1.bam, example2.bam to SJ files (results/example{1,2}.sj)
new-majiq sj bam/example1.bam results/annotated.splicegraph results/example1.sj
new-majiq sj bam/example2.bam results/annotated.splicegraph results/example2.sj

# build to make results/combined.splicegraph
new-majiq build results/annotated.splicegraph results/combined.splicegraph --sjs results/example1.sj results/example2.sj

# get psi coverage
# TODO
```
