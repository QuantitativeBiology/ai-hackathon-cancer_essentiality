"""Model keys of the G3 nonlinear track (registered in scripts/run.py next to essml.models.MODELS).

Kept separate from ``essml.models.__init__`` so the kernel track and this track edit disjoint files.
"""
from __future__ import annotations

from .factorised import FactorisedM4, NeuralPC, RidgePC
from .lgbm import LGBMResidualM3
from .published import ElasticNetR2, TopKCorrR4

MODELS = {"M4pre_MLP": NeuralPC, "M4pre_RidgePC": RidgePC, "M4": FactorisedM4,
          "M3": LGBMResidualM3, "R2": ElasticNetR2, "R4": TopKCorrR4}
