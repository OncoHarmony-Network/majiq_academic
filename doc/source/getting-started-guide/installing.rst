.. _installing:

Installation
============

Required dependencies
---------------------

- Python (3.8 or later)
- setuptools (45 or later)
- HTSLib_ (1.10 or later)
- C++ compiler with C++17 support (e.g. gcc >= 7)
- moccasin >= 0.26 (currently branch new_moccasin_)

.. _HTSLib: https://github.com/samtools/htslib/releases
.. _new_moccasin: https://bitbucket.org/biociphers/moccasin/src/new_moccasin


HTSLib installation
~~~~~~~~~~~~~~~~~~~

MAJIQ requires HTSLib_ (>= 1.10) to be installed.
This may already available in a shared directory of your environment,
otherwise, you can install HTSLib from source.
For example, to install htslib-1.13 to ``~/install/htslib-1.13``, you could run:

.. code-block:: bash

   # download htslib 1.13 archive to current directory
   curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
   tar -xf htslib-1.13.tar.bz2  # extract archive
   cd htslib-1.13  # change into htslib source directory
   # configure, make, and install htslib to ~/install/htslib-1.13
   ./configure --prefix=$HOME/install/htslib-1.13
   make
   make install


See ``INSTALL`` from the HTSLib_ sources for detailed instructions.


Pip installation
~~~~~~~~~~~~~~~~

Install MAJIQ using pip.
MAJIQ setup needs to know where HTSLib is installed.
By default, it assumes that the library and include directories are in the Unix
default locations: ``/usr/local/lib``, ``/usr/local/include``.
If this is not the case, please set the environment variables
``HTSLIB_LIBRARY_DIR`` and ``HTSLIB_INCLUDE_DIR`` to the actual locations.
Install a version of moccasin since the rewrite (>= 0.26). This can currently
be done by running
``pip install git+https://bitbucket.org/biociphers/moccasin@new_moccasin``.
Then install this package from base directory of this repository:
``pip install .``.

For example, with previous instructions, from same directory as README:

.. code-block:: bash

   # change to where library/include directories are installed
   export HTSLIB_LIBRARY_DIR=$HOME/install/htslib-1.13/lib
   export HTSLIB_INCLUDE_DIR=$HOME/install/htslib-1.13/include
   # install moccasin
   pip install git+https://bitbucket.org/biociphers/moccasin@new_moccasin
   # NOTE: from same directory as this README
   pip install .  # install both majiq and voila
