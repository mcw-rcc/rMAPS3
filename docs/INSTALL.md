# Installation

## 1. Python Environment (Recommended)

Using a virtual environment is recommended but not required.

```bash
python -m venv .venv
```

Activate:

- Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

- Linux/macOS:

```bash
source .venv/bin/activate
```

## 2. Install Dependencies

```bash
pip install -r requirements.txt
```

## 3. Prepare Genome FASTA Data

Layout must be:

```text
genomedata/
  <build>/<build>.fa
```

Examples:

```text
genomedata/hg19/hg19.fa
genomedata/hg38/hg38.fa
genomedata/mm10/mm10.fa
```

Download and verify genomes using project scripts:

- Windows: `scripts/fetch_genomes.ps1`
- Linux/macOS: `scripts/fetch_genomes.sh`
- Manifest (source URL + SHA256): `scripts/genomes.manifest.tsv`

## 4. Verify Installation

Quick CLI smoke check:

```bash
python tests/smoke_cli.py
```

Optional command checks:

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## 5. Optional: Start Local Web UI

Run local web UI:

```bash
python run_web.py
```

Then open `http://127.0.0.1:5000`.

## 6. Next Docs

- CLI examples and command reference: `docs/CLI_USAGE.md`
- Web UI details: `webui/README.md`
- Test scripts and matrix: `tests/README.md`
