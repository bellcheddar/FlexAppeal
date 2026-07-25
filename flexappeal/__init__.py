"""FlexAppeal -- prepare OpenMM molecular dynamics runs, analyse what comes back.

The hosted app never simulates. It turns the OpenMM option surface into a guided
form, emits a self-contained bundle you run on your own machine, and then reads
the single results file that run produces.
"""

from .options import FLEXAPPEAL_VERSION, OPENMM_VERSION

__version__ = FLEXAPPEAL_VERSION
__all__ = ["FLEXAPPEAL_VERSION", "OPENMM_VERSION", "__version__"]
