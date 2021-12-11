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

from .core.Contigs import Contigs
from .core.DeltaPsi import DeltaPsi
from .core.DPsiPrior import DPsiPrior
from .core.Events import Events, UniqueEventsMasks
from .core.EventsCoverage import EventsCoverage
from .core.ExonConnections import ExonConnections
from .core.Exons import Exons
from .core.GeneIntrons import GeneIntrons
from .core.GeneJunctions import GeneJunctions
from .core.Genes import Genes
from .core.GroupIntronsGenerator import GroupIntronsGenerator
from .core.Heterogen import Heterogen
from .core.PassedJunctions import GroupJunctionsGenerator, PassedJunctionsGenerator
from .core.PMFSummaries import PMFSummaries

# PsiCoverage and summaries for case/control outliers
from .core.PsiControlsSummary import PsiControlsSummary
from .core.PsiCoverage import PsiCoverage
from .core.PsiOutliers import PsiOutliers
from .core.SimplifierGroup import SimplifierGroup

# how we measure junctions and introns from BAM files
from .core.SJExperiment import SJExperiment
from .core.SJIntrons import SJIntrons
from .core.SJIntronsBins import SJIntronsBins
from .core.SJJunctions import SJJunctions
from .core.SJJunctionsBins import SJJunctionsBins

# model of genes, exons connected by introns/junctions
from .core.SpliceGraph import SpliceGraph

# coverage on splicegraph or events
from .core.SpliceGraphReads import SpliceGraphReads


def rng_seed(seed: int) -> None:
    """Set seed for random number generator pools

    Set seed for random number generator pools. There are two separate rng
    pools in new-majiq, one for nm.beta_mixture, and another for nm.internals.
    This is a helper function that sets the seed for both pools. Note that
    this means that random number generators from these modules could be
    correlated if used in the same session (we offer this convenience function
    because we think that generating LSVCoverage and sampling from beta
    mixtures in the same session is a rare use case).

    Parameters
    ----------
    seed: int
        Integer value used to seed the random number generators
    """
    _internals_rng_seed(seed)
    beta_mixture.rng_seed(seed)
    return


def rng_resize(n: int) -> None:
    """Resize rng pools to allow n simultaneous threads

    Resize rng pools to allow n simultaneous threads. There are two separate
    rng pools in new-majiq, one for nm.beta_mixture, and another for
    nm.internals. This is a helper function that resizes both pools.

    Parameters
    ----------
    n: int
        The number of random number generators that can be simultaneously
        acquired (should be equal to the number of threads performing random
        number generation)
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
    "Contigs",
    "Events",
    "UniqueEventsMasks",
    "EventsCoverage",
    "ExonConnections",
    "Exons",
    "GeneIntrons",
    "GeneJunctions",
    "Genes",
    "GroupIntronsGenerator",
    "GroupJunctionsGenerator",
    "PassedJunctionsGenerator",
    "SimplifierGroup",
    "PsiCoverage",
    "DPsiPrior",
    "DeltaPsi",
    "Heterogen",
    "PassedJunctions",
    "PMFSummaries",
    "PsiControlsSummary",
    "PsiOutliers",
    "SJExperiment",
    "SJIntrons",
    "SJIntronsBins",
    "SJJunctions",
    "SJJunctionsBins",
    "SpliceGraph",
    "SpliceGraphReads",
]
