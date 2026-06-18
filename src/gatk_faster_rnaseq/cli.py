from __future__ import annotations

from collections.abc import Sequence

from .pipeline.runner import main as run_pipeline_main


def main(argv: Sequence[str] | None = None) -> int:
    return run_pipeline_main(argv)

