from new_majiq.internals import ExperimentStrandness, ExperimentThresholds
from new_majiq.internals import rng_resize as _internals_rng_resize
from new_majiq.internals import rng_seed as _internals_rng_seed

try:
    from new_majiq._version import version as __version__
except (ModuleNotFoundError, ImportError):
    __version__ = "3.0.0unknown"
except Exception:
    raise

import new_majiq.beta_mixture as beta_mixture
import new_majiq.constants as constants
import new_majiq.stats as stats

from .DeltaPsi import DeltaPsi
from .DPsiPrior import DPsiPrior
from .Events import Events
from .EventsCoverage import EventsCoverage
from .ExonConnections import ExonConnections
from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .Genes import Genes
from .Heterogen import Heterogen

# summary of PsiCoverage in population
from .PsiControlsSummary import PsiControlsSummary

# coverage in terms of PSI/total coverage on event connections
from .PsiCoverage import PsiCoverage

# how we measure junctions and introns from BAM files
from .SJExperiment import SJExperiment
from .SJIntrons import SJIntrons
from .SJIntronsBins import SJIntronsBins
from .SJJunctions import SJJunctions
from .SJJunctionsBins import SJJunctionsBins

# model of genes, exons connected by introns/junctions
from .SpliceGraph import SpliceGraph

# coverage on splicegraph or events
from .SpliceGraphReads import SpliceGraphReads


def rng_seed(seed: int) -> None:
    """Set seed for random number generator pools

    Set seed for random number generator pools. There are two separate rng
    pools in new-majiq, one for nm.beta_mixture, and another for nm.internals.
    This is a helper function that sets the seed for both pools. Note that
    this means that random number generators from these modules could be
    correlated if used in the same session (we offer this convenience function
    because we think that generating LSVCoverage and sampling from beta
    mixtures in the same session is a rare use case).
    """
    _internals_rng_seed(seed)
    beta_mixture.rng_seed(seed)
    return


def rng_resize(n: int) -> None:
    """Resize rng pools to allow n simultaneous threads

    Resize rng pools to allow n simultaneous threads. There are two separate
    rng pools in new-majiq, one for nm.beta_mixture, and another for
    nm.internals. This is a helper function that resizes both pools.
    """
    _internals_rng_resize(n)
    beta_mixture.rng_resize(n)
    return


__all__ = [
    "beta_mixture",
    "stats",
    "ExperimentStrandness",
    "ExperimentThresholds",
    "rng_resize",
    "rng_seed",
    "__version__",
    "constants",
    "Events",
    "EventsCoverage",
    "ExonConnections",
    "Exons",
    "GeneIntrons",
    "GeneJunctions",
    "Genes",
    "PsiCoverage",
    "DPsiPrior",
    "DeltaPsi",
    "Heterogen",
    "PsiControlsSummary",
    "SJExperiment",
    "SJIntrons",
    "SJIntronsBins",
    "SJJunctions",
    "SJJunctionsBins",
    "SpliceGraph",
    "SpliceGraphReads",
]
