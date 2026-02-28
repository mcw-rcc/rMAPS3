# Installation

## 1. Create Environment

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

## 4. Verify CLI Wiring

```bash
python tests/smoke_cli.py
```
