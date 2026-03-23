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
python -m pip install -r requirements.txt
```

Optional system dependencies (not installed via `pip`):

- MiKTeX or TeX Live: improves PyX text rendering in motif-map outputs
- Ghostscript: enables native PNG export path

## 3. Prepare Genome FASTA Data

Genomes are **not** included in the repository (files are ~3GB each). Download them using the provided scripts:

**Windows (PowerShell):**

```powershell
# Download hg19 and hg38
.\scripts\fetch_genomes.ps1 -Genomes hg19,hg38

# Or all available genomes
.\scripts\fetch_genomes.ps1 -Genomes all

# Dry run to see what would download
.\scripts\fetch_genomes.ps1 -Genomes hg19 -DryRun
```

**Linux/macOS:**

```bash
# Download hg19 and hg38
./scripts/fetch_genomes.sh --genomes hg19,hg38

# Or all available genomes
./scripts/fetch_genomes.sh --genomes all

# Dry run to see what would download
./scripts/fetch_genomes.sh --genomes hg19 --dry-run
```

The script:
- Downloads from UCSC/Ensembl sources
- Verifies SHA256 checksums automatically
- Decompresses `.fa.gz` → `.fa`
- Creates `genomedata/<build>/<build>.fa` structure

**Available genomes:** hg19, hg38, mm10, dm3, dm6, and others (see `scripts/genomes.manifest.tsv`)

## 4. Verify Installation

Quick CLI smoke check:

```bash
python tests/smoke_cli.py
```

Recommended pre-deployment validation:

```bash
python tests/test_clip.py
python tests/test_motif.py --fasta-root genomedata --genome hg19
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

- CLI examples and command reference: [`CLI_USAGE.md`](CLI_USAGE.md)
- Web UI details: [`../webui/README.md`](../webui/README.md)
- Test scripts and matrix: [`../tests/README.md`](../tests/README.md)
