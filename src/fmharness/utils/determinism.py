# pyright: reportMissingImports=none, reportUnknownMemberType=none
"""Deterministic-execution helpers.

Call ``fix_seeds(seed)`` from every CLI entrypoint before sampling, model
loading, or inference. Pairs with ``EnvironmentSnapshot.cuda_deterministic``.
"""

from __future__ import annotations

import os
import random

import numpy as np


def fix_seeds(seed: int) -> None:
    """Set deterministic seeds and CUDA workspace config across the stack.

    ``CUBLAS_WORKSPACE_CONFIG`` must be set before torch initializes CUDA
    state, so callers must invoke this before importing or using any torch
    CUDA functionality.
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch
    except ImportError:
        return

    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
