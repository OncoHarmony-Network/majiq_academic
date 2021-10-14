from typing import Dict

import numpy

approximation: numpy.ufunc
cdf: numpy.ufunc
dpsi_discrete: numpy.ufunc
moments: numpy.ufunc
pdf: numpy.ufunc
pmf: numpy.ufunc
quantile: numpy.ufunc
sample: numpy.ufunc
ttest_sample: numpy.ufunc
stats_sample: numpy.ufunc

def rng_resize(n: int) -> None: ...
def rng_seed(seed: int) -> None: ...
def stats_available() -> Dict[str, int]: ...
