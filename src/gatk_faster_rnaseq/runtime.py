from __future__ import annotations

import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class Step:
    name: str
    command: list[str]
    outputs: tuple[Path, ...]
    env: dict[str, str] | None = None


def quote_cmd(command: list[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def append_timing(path: Path, step: str, status: str, seconds: float, command: str, message: str) -> None:
    new_file = not path.exists()
    with path.open("a") as handle:
        if new_file:
            handle.write("step\tstatus\tseconds\tminutes\tcommand\tmessage\n")
        handle.write(f"{step}\t{status}\t{seconds:.3f}\t{seconds / 60:.3f}\t{command}\t{message}\n")


def run_step(step: Step, log_dir: Path, timings_path: Path, force: bool) -> None:
    if step.outputs and not force and all(path.exists() for path in step.outputs):
        append_timing(timings_path, step.name, "skipped", 0.0, quote_cmd(step.command), "")
        return

    log_path = log_dir / f"{step.name}.log"
    start = time.perf_counter()
    status = "ok"
    message = ""
    with log_path.open("w") as log:
        _write_command(log, step.command)
        try:
            subprocess.run(
                step.command,
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=_merged_env(step.env),
            )
        except subprocess.CalledProcessError as exc:
            status = "failed"
            message = f"exit_code={exc.returncode}; log={log_path}"
            raise
        finally:
            elapsed = time.perf_counter() - start
            append_timing(timings_path, step.name, status, elapsed, quote_cmd(step.command), message)


def run_parallel_steps(steps: list[Step], log_dir: Path, timings_path: Path, force: bool, max_parallel: int) -> None:
    if max_parallel < 1:
        raise ValueError("max_parallel must be at least 1")

    for batch_start in range(0, len(steps), max_parallel):
        running: list[tuple[Step, subprocess.Popen[bytes], TextIO, Path, float]] = []
        first_error: subprocess.CalledProcessError | None = None

        for step in steps[batch_start : batch_start + max_parallel]:
            if step.outputs and not force and all(path.exists() for path in step.outputs):
                append_timing(timings_path, step.name, "skipped", 0.0, quote_cmd(step.command), "")
                continue

            log_path = log_dir / f"{step.name}.log"
            log = log_path.open("w")
            _write_command(log, step.command)
            start = time.perf_counter()
            process = subprocess.Popen(
                step.command,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=_merged_env(step.env),
            )
            running.append((step, process, log, log_path, start))

        for step, process, log, log_path, start in running:
            status = "ok"
            message = ""
            return_code = process.wait()
            log.close()
            if return_code != 0:
                status = "failed"
                message = f"exit_code={return_code}; log={log_path}"
                if first_error is None:
                    first_error = subprocess.CalledProcessError(return_code, step.command)
            elapsed = time.perf_counter() - start
            append_timing(timings_path, step.name, status, elapsed, quote_cmd(step.command), message)

        if first_error is not None:
            raise first_error


def _write_command(log: TextIO, command: list[str]) -> None:
    log.write(f"$ {quote_cmd(command)}\n\n")
    log.flush()


def _merged_env(overrides: dict[str, str] | None) -> dict[str, str] | None:
    if overrides is None:
        return None
    env = os.environ.copy()
    env.update(overrides)
    return env
