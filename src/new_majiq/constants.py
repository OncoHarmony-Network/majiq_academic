"""
constants.py

Constants for new_majiq front-end

Author: Joseph K Aicher
"""

import string
from typing import Dict, Final, List, Set

from new_majiq.beta_mixture import stats_available as _stats_available
from new_majiq.internals import (
    ExperimentStrandness,
    ExperimentThresholds,
    GFF3FeatureType,
    _default_gff3_types,
)

NC_CONTIGS: Final[str] = "contigs"
NC_GENES: Final[str] = "genes"
NC_EXONS: Final[str] = "exons"
NC_GENEJUNCTIONS: Final[str] = "junctions"
NC_GENEINTRONS: Final[str] = "introns"

NC_SJINTRONSCONTIGS: Final[str] = "sj_introns_contigs"
NC_SJINTRONS: Final[str] = "sj_introns"
NC_SJINTRONSBINS: Final[str] = "sj_introns_bins"
NC_SJJUNCTIONS: Final[str] = "sj_junctions"
NC_SJJUNCTIONSBINS: Final[str] = "sj_junctions_bins"
NC_SJJUNCTIONSCONTIGS: Final[str] = "sj_junctions_contigs"

NC_EVENTS: Final[str] = "events"
NC_EVENTSCOVERAGE: Final[str] = "events_coverage"
NC_PSICOVERAGE: Final[str] = "psi_coverage"

NC_EVENTSQUANTIFIED: Final[str] = "events_quantified"

NC_SGREADS: Final[str] = "sg_reads"
# number of junctions/introns per chunk in sg_reads
# size of computation across N experiments for a chunk is
# N * 8 * chunksize bytes. So N=20000, chunksize = 1<<15 is under 5GB
NC_SGREADS_CHUNKS: Final[int] = 1 << 15


DEFAULT_BUILD_PROCESS_IR: Final[bool] = True
DEFAULT_BUILD_GFF3TYPES: Final[Dict[str, GFF3FeatureType]] = _default_gff3_types()

DEFAULT_BAM_STRANDNESS: Final[ExperimentStrandness] = ExperimentStrandness.NONE
DEFAULT_BAM_NTHREADS: Final[int] = 1
DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS: Final[bool] = False
# how we detect strandness
DEFAULT_BAM_STRAND_MINREADS: Final[int] = 10
DEFAULT_BAM_STRAND_MINJUNCTIONS: Final[int] = 100
DEFAULT_BAM_STRAND_MINDEVIATION: Final[float] = 0.2

DEFAULT_BUILD_MINREADS: Final[int] = 3
DEFAULT_BUILD_MINDENOVO: Final[int] = 5
DEFAULT_BUILD_MINPOS: Final[int] = 2
DEFAULT_BUILD_MAX_PCTBINS: Final[float] = 0.6
DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY: Final[float] = 0.5
DEFAULT_BUILD_MATCH_INTRON_PROBABILITY: Final[float] = 0.95
DEFAULT_BUILD_EXP_THRESHOLDS: Final[ExperimentThresholds] = ExperimentThresholds(
    minreads=DEFAULT_BUILD_MINREADS,
    mindenovo=DEFAULT_BUILD_MINDENOVO,
    minpos=DEFAULT_BUILD_MINPOS,
    max_pctbins=DEFAULT_BUILD_MAX_PCTBINS,
    junction_acceptance_probability=DEFAULT_BUILD_MATCH_JUNCTION_PROBABILITY,
    intron_acceptance_probability=DEFAULT_BUILD_MATCH_INTRON_PROBABILITY,
)
DEFAULT_BUILD_MINEXPERIMENTS: Final[float] = 0.5
DEFAULT_BUILD_DENOVO_JUNCTIONS: Final[bool] = True
DEFAULT_BUILD_DENOVO_IR: Final[bool] = True
DEFAULT_BUILD_KEEP_ANNOTATED_IR: Final[bool] = True
DEFAULT_BUILD_DENOVO_SIMPLIFIED: Final[bool] = False

DEFAULT_SIMPLIFIER_MINEXPERIMENTS: Final[float] = 1.0
DEFAULT_SIMPLIFIER_MINPSI: Final[float] = 0.01
DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED: Final[float] = 0.0
DEFAULT_SIMPLIFIER_MINREADS_DENOVO: Final[float] = 0.0
DEFAULT_SIMPLIFIER_MINREADS_INTRON: Final[float] = 0.0

DEFAULT_COVERAGE_NUM_BOOTSTRAPS: Final[int] = 30
DEFAULT_COVERAGE_STACK_PVALUE: Final[float] = 1e-7

DEFAULT_QUANTIFY_NTHREADS: Final[int] = 1
DEFAULT_QUANTIFY_MINREADS: Final[float] = 10.0
DEFAULT_QUANTIFY_MINBINS: Final[float] = 3.0
DEFAULT_QUANTIFY_MINEXPERIMENTS: Final[float] = 0.5
DEFAULT_QUANTIFY_PSIBINS: Final[int] = 40

DEFAULT_DPSI_PRIOR_MINREADS: Final[float] = 3.0 * DEFAULT_QUANTIFY_MINREADS
DEFAULT_DPSI_PRIOR_MINLSV: Final[int] = 100
DEFAULT_DPSI_PRIOR_MAXITER: Final[int] = 1
DEFAULT_DPSI_PRIOR_A: Final[List[float]] = [1.0, 75.0, 1000.0]
DEFAULT_DPSI_PRIOR_PMIX: Final[List[float]] = [0.2, 0.5, 0.3]
DEFAULT_DPSI_CHANGING_THRESHOLD: Final[float] = 0.20
DEFAULT_DPSI_NONCHANGING_THRESHOLD: Final[float] = 0.05

STATS_AVAILABLE: Final[Dict[str, int]] = _stats_available()

ALLOWED_GROUP_NAME_CHARS: Final[Set[str]] = {*string.ascii_letters, *string.digits, "_"}

DEFAULT_MOCCASIN_RUV_MAX_EVENTS: Final[int] = 10000
DEFAULT_MOCCASIN_RUV_MAX_FACTORS: Final[int] = 0

PSICOV_APPROX: Final[Set[str]] = {"approximation", "both"}
PSICOV_BOOTSTRAP: Final[Set[str]] = {"bootstrap", "both"}
PSICOV_POSTERIORS: Final[Set[str]] = PSICOV_APPROX | PSICOV_BOOTSTRAP
DEFAULT_PSICOV_POSTERIOR: Final[str] = "approximation"
assert DEFAULT_PSICOV_POSTERIOR in PSICOV_POSTERIORS

DPSI_SMOOTH: Final[Set[str]] = {"smooth", "both"}
DPSI_LEGACY: Final[Set[str]] = {"legacy", "both"}
DPSI_POSTERIORS: Final[Set[str]] = DPSI_SMOOTH | DPSI_LEGACY
DEFAULT_DPSI_POSTERIOR: Final[str] = "smooth"
assert DEFAULT_DPSI_POSTERIOR in DPSI_POSTERIORS
