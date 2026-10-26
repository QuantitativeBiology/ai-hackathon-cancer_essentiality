"""Model keys of the G3 nonlinear track (registered in scripts/run.py next to essml.models.MODELS).

Kept separate from ``essml.models.__init__`` so the kernel track and this track edit disjoint files.
"""
from __future__ import annotations

from .factorised import FactorisedM4, NeuralPC, RidgePC

MODELS = {"M4pre_MLP": NeuralPC, "M4pre_RidgePC": RidgePC, "M4": FactorisedM4}

try:  # optional rungs (added as they are implemented)
    from .lgbm import LGBMResidualM3
    MODELS["M3"] = LGBMResidualM3
except ImportError:  # pragma: no cover
    pass
try:
    from .published import ElasticNetR2, TopKCorrR4
    MODELS["R2"] = ElasticNetR2
    MODELS["R4"] = TopKCorrR4
except ImportError:  # pragma: no cover
    pass
