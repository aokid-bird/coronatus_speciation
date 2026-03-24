#!/usr/bin/env python3

import argparse
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run one TreeMix replicate with timeout and retry handling."
    )
    parser.add_argument("--input", required=True, help="Compressed TreeMix input (.gz)")
    parser.add_argument("--outprefix", required=True, help="TreeMix output prefix")
    parser.add_argument("--edge", required=True, type=int, help="Migration edge count")
    parser.add_argument("--block", required=True, type=int, help="TreeMix -k block size")
    parser.add_argument("--timeout-seconds", required=True, type=int, help="Per-attempt timeout in seconds")
    parser.add_argument("--max-attempts", required=True, type=int, help="Maximum retry attempts")
    parser.add_argument("--root-opt", default="", help="Optional TreeMix root argument string")
    parser.add_argument("--se-flag", default="", help="Optional TreeMix standard error flag")
    parser.add_argument("--bootstrap", action="store_true", help="Pass -bootstrap to TreeMix")
    parser.add_argument("--log", required=True, help="Log file path")
    return parser.parse_args()


def cleanup_outputs(prefix: str):
    parent = Path(prefix).parent
    stem = Path(prefix).name
    for path in parent.glob(f"{stem}*"):
        if path.is_file() or path.is_symlink():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def main():
    args = parse_args()
    timeout_seconds = max(1, args.timeout_seconds)
    max_attempts = max(1, args.max_attempts)
    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "treemix",
        "-i", args.input,
        "-m", str(args.edge),
        "-o", str(outprefix),
        "-k", str(args.block),
    ]
    if args.bootstrap:
        cmd.append("-bootstrap")
    if args.root_opt.strip():
        cmd.extend(args.root_opt.strip().split())
    if args.se_flag.strip():
        cmd.extend(args.se_flag.strip().split())

    llik_path = Path(f"{outprefix}.llik")

    for attempt in range(1, max_attempts + 1):
        cleanup_outputs(str(outprefix))
        with log_path.open("a", encoding="utf-8") as log_handle:
            log_handle.write(
                f"[treemix] attempt={attempt}/{max_attempts} edge={args.edge} "
                f"timeout_seconds={timeout_seconds}\n"
            )
            log_handle.write("[treemix] command: " + " ".join(cmd) + "\n")
            log_handle.flush()
            proc = subprocess.Popen(
                cmd,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                text=True,
            )
            try:
                rc = proc.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                log_handle.write(
                    f"[treemix] timeout after {timeout_seconds}s; killing process group.\n"
                )
                log_handle.flush()
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    proc.wait()
                continue

            if rc == 0 and llik_path.exists() and llik_path.stat().st_size > 0:
                log_handle.write("[treemix] completed successfully.\n")
                return 0

            log_handle.write(
                f"[treemix] failed with exit code {rc}; expected output {llik_path} not ready.\n"
            )

    cleanup_outputs(str(outprefix))
    print(
        f"TreeMix failed for edge={args.edge} after {max_attempts} attempts.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
