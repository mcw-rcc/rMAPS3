# rMAPS Refactored

A Python 3 refactor of the original rMAPS motif-mapping pipeline.

## What This Project Provides

- Python 3 support across legacy motif-map event engines.
- `pyfaidx`-based genome access (no `pygr` runtime dependency).
- Unified CLI (`cli.py`) for:
  - motif-map generation (`se`, `a3ss`, `a5ss`, `ri`, `mxe`)
  - MISO-to-rMATS conversion
  - exon-set generation
- Local Web UI (`run_web.py` / `webui/`) for browser-driven runs.

## Requirements

- Python 3.10+ recommended
- OS: Linux/macOS/Windows
- Perl available on `PATH` for MISO conversion (`miso2rMATS*.pl`)
- Flask for Web UI (`run_web.py`)
- Optional Ghostscript for native PNG export (not required with Python fallback)

Install dependencies:

```bash
pip install -r requirements.txt
```

## Documentation

- Installation: [`docs/INSTALL.md`](docs/INSTALL.md)
- CLI reference: [`docs/CLI_USAGE.md`](docs/CLI_USAGE.md)
- Web UI usage and API: [`webui/README.md`](webui/README.md)
- Test scripts and matrix: [`tests/README.md`](tests/README.md)

## Genome Data Layout

This repo expects a FASTA root directory passed as `--fasta-root` (or `--fastaRoot`) with:

```text
genomedata/
  hg19/hg19.fa
  hg38/hg38.fa
  mm10/mm10.fa
  dm3/dm3.fa
  ...
```

Important:
- Build name and FASTA name must match (`<build>/<build>.fa`).
- If you run with `--genome hg19`, the loader uses `genomedata/hg19/hg19.fa`.

For reproducible genome setup (download scripts + SHA256 manifest), see [`docs/INSTALL.md`](docs/INSTALL.md).

## Project Structure

- `cli.py`: primary CLI entrypoint
- `run_web.py`: local web launcher entrypoint
- `rmaps_core/`: core importable modules
- `legacy/`: event-specific motif scripts (`motifMap*_MP.py`)
- `bin/`: converter/util scripts
- `scripts/`: genome download and setup helpers
- `data/`: motif tables and reference motif inputs
- `testData/`: bundled sample inputs for quick tests
- `genomedata/`: local genome FASTA root (`<build>/<build>.fa`)
- `docs/`: installation and CLI docs
- `tests/`: smoke and integration test runners
- `webui/`: Flask app, static assets, and templates

## CLI Usage

Top-level help:

```bash
python cli.py --help
```

Motif-map help:

```bash
python cli.py motif-map --help
python cli.py motif-map se --help
```

Converters:

```bash
python cli.py convert --help
python cli.py convert miso --help
```

Exon sets:

```bash
python cli.py exon-sets --help
python cli.py exon-sets se --help
```

Detailed command examples are in [`docs/CLI_USAGE.md`](docs/CLI_USAGE.md).

## Web UI (Local)

Run the local web server:

```bash
python run_web.py
```

Open:

- `http://127.0.0.1:5000`

See full Web UI documentation in [`webui/README.md`](webui/README.md).

## Outputs

Each run writes to the given `--output` directory:

- `exon/`: `up.coord.txt`, `dn.coord.txt`, `bg.coord.txt`
- `fasta/`: per-region FASTA windows
- `maps/`: motif map outputs (`.pdf`, `.png` where available)
- `temp/`: intermediate files
- `log.motifMap.txt`: run log

## Testing

Detailed test documentation and commands are in [`tests/README.md`](tests/README.md).

Quick smoke:

```bash
python tests/smoke_cli.py
```

Full event matrix:

```bash
bash tests/run_all_events.sh
```

Note: shell test runners under `tests/` require a Bash-compatible environment (Linux/macOS, WSL, or Git Bash on Windows).

## Operational Notes

- Use the CLI (`cli.py`) as the public interface.
- For browser-based runs, use `run_web.py`.
- Core modules are under `rmaps_core/`.
- Scripts under `legacy/` are internal implementation details; use the CLI for production runs.
- `--fasta-root` is the canonical argument; `--fastaRoot` is accepted for compatibility.

## Troubleshooting

- `FastaNotFoundError`:
  - Verify `--genome` and `--fasta-root` match expected layout (`<build>/<build>.fa`).
- MISO conversion fails:
  - Ensure Perl is installed and available in `PATH`.
  - Verify `bin/miso2rMATS.*.pl` scripts exist.
- PNG files missing:
  - Native Ghostscript export may be unavailable.
  - Install optional rendering dependencies from `requirements.txt`.
