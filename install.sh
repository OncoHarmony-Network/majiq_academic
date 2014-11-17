#!/bin/bash

if [ "$#" -ne 2 ]
then
  echo " This scripts adapts the compiled libraries to your computer. It has one parameter, "
  echo " the absolute path where the libraries and binaries will be located "
  echo "Usage: ./install.sh <path to binaries destinations> "
  exit 1
fi

chmod u+w $1/*dylib

#install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "\$1/libgfortran.3.dylib" $1/_ufuncs.so 
#install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgcc_s.1.dylib" "$1/libgcc_s.1.dylib" $1/_ufuncs.so 
#install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/_ufuncs.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.special._ufuncs.so
install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgcc_s.1.dylib" "$1/libgcc_s.1.dylib" $1/scipy.special._ufuncs.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.special._ufuncs.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/libgfortran.3.dylib
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/libgfortran.3.dylib
install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/libgfortran.3.dylib
install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgcc_s.1.dylib" "$1/libgcc_s.1.dylib" $1/libgfortran.3.dylib 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/libgfortran.3.dylib 
install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgcc_s.1.dylib" "$1/libgcc_s.1.dylib" $1/libquadmath.0.dylib 
install_name_tool -change "/usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgcc_s.1.dylib" "$1/libgcc_s.1.dylib" $1/libquadmath.0.dylib 
install_name_tool -change "/usr/local/lib/libfreetype.6.dylib" "$1/libfreetype.6.dylib" $1/matplotlib.ft2font.so 
install_name_tool -change "/usr/local/lib/libpng16.16.dylib" "$1/libpng16.16.dylib" $1/matplotlib.ft2font.so 
install_name_tool -change "/usr/local/lib/libpng16.16.dylib" "$1/libpng16.16.dylib" $1/libfreetype.6.dylib
install_name_tool -change "/usr/local/lib/libpng16.16.dylib" "$1/libpng16.16.dylib" $1/matplotlib._png.so
install_name_tool -change "/usr/local/lib/libfreetype.6.dylib" "$1/libfreetype.6.dylib" $1/matplotlib.backends._backend_agg.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.special.specfun.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.special.specfun.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.linalg._fblas.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.linalg._fblas.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.linalg._flapack.so 
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.linalg._flapack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.linalg.calc_lwork.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.linalg.calc_lwork.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.stats.futil.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.stats.futil.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize.minpack2.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize.minpack2.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize._lbfgsb.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize._lbfgsb.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize._cobyla.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize._cobyla.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize._slsqp.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize._slsqp.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize._minpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize._minpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.sparse.linalg.isolve._iterative.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.sparse.linalg.isolve._iterative.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.sparse.linalg.eigen.arpack._arpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.sparse.linalg.eigen.arpack._arpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.optimize._nnls.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.optimize._nnls.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.integrate._odepack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.integrate._odepack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.integrate._quadpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.integrate._quadpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.integrate.vode.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.integrate.vode.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.integrate._dop.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.integrate._dop.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.integrate.lsoda.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.integrate.lsoda.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.stats.statlib.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.stats.statlib.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.stats.mvn.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.stats.mvn.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.interpolate._fitpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.interpolate._fitpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libgfortran.3.dylib" "$1/libgfortran.3.dylib" $1/scipy.interpolate.dfitpack.so
install_name_tool -change "/usr/local/lib/gcc/x86_64-apple-darwin12.5.0/4.9.1/libquadmath.0.dylib" "$1/libquadmath.0.dylib" $1/scipy.interpolate.dfitpack.so
