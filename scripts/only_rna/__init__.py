"""OnlyRNA config package (minimal foundation).

Expose the core config symbols for tests and downstream imports.
"""

from .models import RunConfig, QcThresholds, PlottingConfig  # noqa: F401
from .config import load_default_config, merge_cli_overrides  # noqa: F401

__all__ = [
    "RunConfig",
    "QcThresholds",
    "PlottingConfig",
    "load_default_config",
    "merge_cli_overrides",
]
