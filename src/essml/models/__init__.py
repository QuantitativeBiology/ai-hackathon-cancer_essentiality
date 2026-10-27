"""Model rungs beyond the baselines (DESIGN.md §5). ``MODELS`` extends ``essml.baselines.MODELS``."""
from __future__ import annotations

from .kernel_ridge import KernelRidgeM1, KernelRidgeR1, KernelRidgeT
from .m2_local import LocalStageM2
from .m5_stack import StackM5

MODELS = {"M1": KernelRidgeM1, "M1b": KernelRidgeM1, "R1": KernelRidgeR1, "T": KernelRidgeT,
          "M2": LocalStageM2, "M5": StackM5}
