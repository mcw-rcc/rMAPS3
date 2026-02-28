"""
Core dispatcher for CLIP-seq RNA map analysis.

This module provides the dispatcher logic for running CLIP-seq analysis
across all 5 alternative splicing event types (SE, A3SS, A5SS, RI, MXE).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import sys
import subprocess
import os


PYTHON = sys.executable
REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ClipEventSpec:
    """
    Specification for each CLIP-seq event type.

    Attributes:
        name: Full event name (e.g., "SE", "A5SS")
        script: Path to the legacy clip processing script
        description: Human-readable description of the event type
    """

    name: str
    script: str
    description: str


CLIP_EVENT_SPECS: Dict[str, ClipEventSpec] = {
    "se": ClipEventSpec(
        name="SE",
        script="legacy/clipSeqSE.py",
        description="Skipped Exon (Exon Skipping)",
    ),
    "a3ss": ClipEventSpec(
        name="A3SS",
        script="legacy/clipSeqA3SS.py",
        description="Alternative 3' Splice Site",
    ),
    "a5ss": ClipEventSpec(
        name="A5SS",
        script="legacy/clipSeqA5SS.py",
        description="Alternative 5' Splice Site",
    ),
    "ri": ClipEventSpec(
        name="RI",
        script="legacy/clipSeqRI.py",
        description="Retained Intron",
    ),
    "mxe": ClipEventSpec(
        name="MXE",
        script="legacy/clipSeqMXE.py",
        description="Mutually Exclusive Exons",
    ),
}


def clip_event_script(event: str) -> Path:
    """
    Return the path to the CLIP-seq processing script for the given event.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)

    Returns:
        Path to the clip processing script

    Raises:
        ValueError: If event type is not supported
    """
    key = event.lower()
    if key not in CLIP_EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return REPO_ROOT / CLIP_EVENT_SPECS[key].script


def run_subprocess(cmd: list[str]) -> int:
    """
    Shared helper for launching child processes from this project.

    Args:
        cmd: Command list to execute

    Returns:
        Return code from the subprocess
    """
    env = dict(os.environ)
    repo_root = str(REPO_ROOT)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = repo_root if not existing else f"{repo_root}{os.pathsep}{existing}"
    result = subprocess.run(cmd, cwd=REPO_ROOT, env=env)
    return result.returncode


def run_clip_map(
    event: str,
    peak: Path,
    output: Path,
    rmats: str,
    miso: str,
    up: str,
    down: str,
    background: str,
    label: str,
    intron: int,
    exon: int,
    window: int,
    step: int,
    sig_fdr: float,
    sig_delta_psi: float,
    separate: bool,
) -> int:
    """
    Run CLIP-seq RNA map analysis for the specified event type.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)
        peak: Path to CLIP-seq peak file
        output: Output directory path
        rmats: rMATS file path or "NA"
        miso: MISO file path or "NA"
        up: Upregulated exons coordinate file or "NA"
        down: Downregulated exons coordinate file or "NA"
        background: Background exons coordinate file or "NA"
        label: RNA-binding protein label
        intron: Intron length to examine (default: 250)
        exon: Exon length to examine (default: 50)
        window: Window size for analysis (default: 10)
        step: Step size for sliding window (default: 1)
        sig_fdr: FDR cutoff for significant events
        sig_delta_psi: Delta PSI cutoff for significant events
        separate: Whether to create separate p-value plots

    Returns:
        Return code from the subprocess
    """
    script = clip_event_script(event)

    cmd = [
        PYTHON,
        str(script),
        "--peak", str(peak),
        "--output", str(output),
        "--rMATS", rmats,
        "--miso", miso,
        "--up", up,
        "--down", down,
        "--background", background,
        "--label", label,
        "--intron", str(intron),
        "--exon", str(exon),
        "--window", str(window),
        "--step", str(step),
        "--sigFDR", str(sig_fdr),
        "--sigDeltaPSI", str(sig_delta_psi),
    ]

    if separate:
        cmd.append("--separate")

    return run_subprocess(cmd)


def get_event_description(event: str) -> str:
    """
    Get the human-readable description for an event type.

    Args:
        event: Event type (se, a3ss, a5ss, ri, mxe)

    Returns:
        Description string

    Raises:
        ValueError: If event type is not supported
    """
    key = event.lower()
    if key not in CLIP_EVENT_SPECS:
        raise ValueError(f"Unsupported event type: {event}")
    return CLIP_EVENT_SPECS[key].description


def list_supported_events() -> list[str]:
    """
    Get list of all supported event types.

    Returns:
        List of event type keys (lowercase)
    """
    return list(CLIP_EVENT_SPECS.keys())
