# Testing Guide

Two unified scripts cover routine testing:

- `tests/test_clip.py`
- `tests/test_motif.py`

## Quick Start

```bash
# CLI smoke checks
python tests/smoke_cli.py

# CLIP map test suite
python tests/test_clip.py

# Motif map test suite (requires FASTA data)
python tests/test_motif.py --fasta-root genomedata --genome hg19
```

## CLIP Test Suite

Script: `tests/test_clip.py`

Usage:

```bash
python tests/test_clip.py
python tests/test_clip.py --event se a5ss
python tests/test_clip.py --verbose
python tests/test_clip.py --list
```

Event types:

- `se`
- `a5ss`
- `a3ss`
- `mxe`
- `ri`

## Motif Test Suite

Script: `tests/test_motif.py`

Usage:

```bash
python tests/test_motif.py
python tests/test_motif.py --event se a3ss
python tests/test_motif.py --fasta-root genomedata --genome hg19
python tests/test_motif.py --verbose
python tests/test_motif.py --list
```

Event types:

- `se`
- `a5ss`
- `a3ss`
- `mxe`
- `ri`

## Prerequisites

Common:

- Python 3.10+
- Dependencies from `requirements.txt`
- Run tests from repository root

CLIP tests use:

- `data/test/clip/PIPE-CLIP.Clusters.bed`
- `data/test/clip/ES/se.rMATS.txt`
- `data/test/clip/A3SS/a3ss.rMATS.txt`
- `data/test/clip/A5SS/a5ss.rMATS.txt`
- `data/test/clip/RI/ri.rMATS.txt`
- `data/test/clip/MXE/mxe.rMATS.txt`

Motif tests use:

- `data/testMotifs.txt`
- `data/test/clip/ES/se.rMATS.txt`
- `data/test/clip/A3SS/a3ss.rMATS.txt`
- `data/test/clip/A5SS/a5ss.rMATS.txt`
- `data/test/clip/RI/ri.rMATS.txt`
- `data/test/clip/MXE/mxe.rMATS.txt`
- FASTA files under `<fasta-root>/<genome>/<genome>.fa`

Optional:

- Perl (MISO conversion)
- Ghostscript (PNG generation)

## Legacy Shell Scripts

Legacy shell runners are kept under `tests/legacy/` for backward compatibility.
Use the Python suites above for standard validation.


