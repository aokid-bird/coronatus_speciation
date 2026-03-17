#!/usr/bin/env python3
import shutil
import subprocess
import sys


SUCCESS_STATES = {
    "BOOT_FAIL": "failed",
    "CANCELLED": "failed",
    "COMPLETED": "success",
    "COMPLETING": "running",
    "CONFIGURING": "running",
    "DEADLINE": "failed",
    "FAILED": "failed",
    "NODE_FAIL": "failed",
    "OUT_OF_MEMORY": "failed",
    "PENDING": "running",
    "PREEMPTED": "failed",
    "RUNNING": "running",
    "RESIZING": "running",
    "REQUEUED": "running",
    "SUSPENDED": "running",
    "TIMEOUT": "failed",
}


def _run(cmd):
    try:
        return subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError:
        return None


def _normalize_state(raw_state):
    state = (raw_state or "").strip()
    if not state:
        return None
    state = state.split()[0].split("+")[0]
    return SUCCESS_STATES.get(state, "running")


def _status_from_sacct(jobid):
    if shutil.which("sacct") is None:
        return None
    proc = _run(
        [
            "sacct",
            "-j",
            jobid,
            "--format=State",
            "--noheader",
            "--parsable2",
        ]
    )
    if proc is None or proc.returncode != 0:
        return None

    for line in reversed(proc.stdout.splitlines()):
        state = _normalize_state(line.split("|", 1)[0])
        if state:
            return state
    return None


def _status_from_squeue(jobid):
    if shutil.which("squeue") is None:
        return None
    proc = _run(["squeue", "-h", "-j", jobid, "-o", "%T"])
    if proc is None or proc.returncode != 0:
        return None

    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    if not lines:
        return "success"
    return _normalize_state(lines[0]) or "running"


def main():
    if len(sys.argv) != 2:
        print("running")
        return 0

    jobid = sys.argv[1]
    status = _status_from_sacct(jobid)
    if status is None:
        status = _status_from_squeue(jobid)
    if status is None:
        status = "running"

    print(status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
