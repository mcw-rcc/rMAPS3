# CLI Usage

This document is the complete CLI reference for rMAPS 3.

## Global Help

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py clip-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## Half 1: Motif Map

### Supported Event Types

- `se`: Skipped Exon
- `a3ss`: Alternative 3' Splice Site
- `a5ss`: Alternative 5' Splice Site
- `ri`: Retained Intron
- `mxe`: Mutually Exclusive Exons

### Event Subcommand Help

```bash
python cli.py motif-map se --help
python cli.py motif-map a3ss --help
python cli.py motif-map a5ss --help
python cli.py motif-map ri --help
python cli.py motif-map mxe --help
```

### Shared Required Options

- `--known-motifs`
- `--fasta-root`
- `--genome`
- `--output`

### Input Modes (choose one)

- rMATS mode: `--rMATS <file>` with `--miso NA --up NA --down NA --background NA`
- MISO mode: `--miso <file>` with `--rMATS NA --up NA --down NA --background NA`
- Coordinate replay mode: `--up <file> --down <file> --background <file>` with `--rMATS NA --miso NA`

### Common Optional Options

- `--motifs` (default `NA`)
- `--label` (default `RBP`)
- `--intron` (default `250`)
- `--exon` (default `50`)
- `--window` (default `50`)
- `--step` (default `1`)
- `--sig-fdr` / `--sigFDR` (default `0.05`)
- `--sig-delta-psi` / `--sigDeltaPSI` (default `0.05`)
- `--separate`

### Event Details and Examples

#### `se` (Skipped Exon)
Use when your input file contains cassette-exon skipping events.

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_se \
  --rMATS data/test/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA --up NA --down NA --background NA
```

#### `a3ss` (Alternative 3' Splice Site)
Use for alternative acceptor events.

```bash
python cli.py motif-map a3ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_a3ss \
  --rMATS data/test/clip/A3SS/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `a5ss` (Alternative 5' Splice Site)
Use for alternative donor events.

```bash
python cli.py motif-map a5ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_a5ss \
  --rMATS data/test/clip/A5SS/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `ri` (Retained Intron)
Use for intron-retention events.

```bash
python cli.py motif-map ri \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_ri \
  --rMATS data/test/clip/RI/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

#### `mxe` (Mutually Exclusive Exons)
Use for mutually exclusive exon events.

```bash
python cli.py motif-map mxe \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/motif_mxe \
  --rMATS data/test/clip/MXE/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA
```

## Half 2: CLIP Map

### Supported Event Types

- `se`: Skipped Exon
- `a3ss`: Alternative 3' Splice Site
- `a5ss`: Alternative 5' Splice Site
- `ri`: Retained Intron
- `mxe`: Mutually Exclusive Exons

### Event Subcommand Help

```bash
python cli.py clip-map se --help
python cli.py clip-map a3ss --help
python cli.py clip-map a5ss --help
python cli.py clip-map ri --help
python cli.py clip-map mxe --help
```

### Shared Required Options

- `--peak` / `-p`
- `--output` / `-o`

### Input Modes (choose one)

- rMATS mode: `--rMATS <file>` with `--miso NA --up NA --down NA --background NA`
- MISO mode: `--miso <file>` with `--rMATS NA --up NA --down NA --background NA`
- Coordinate replay mode: `--up <file> --down <file> --background <file>` with `--rMATS NA --miso NA`

### Common Optional Options

- `--label` (default `RBP`)
- `--intron` (default `250`)
- `--exon` (default `50`)
- `--window` (default `10`)
- `--step` (default `1`)
- `--sig-fdr` / `--sigFDR`
- `--sig-delta-psi` / `--sigDeltaPSI`
- `--separate`

Default thresholds by event:

- `se`, `a3ss`, `ri`, `mxe`: `sigFDR=0.05`, `sigDeltaPSI=0.05`
- `a5ss`: `sigFDR=0.005`, `sigDeltaPSI=0.01`

### Event Details and Examples

#### `se` (Skipped Exon)
Use for cassette-exon CLIP enrichment maps.

```bash
python cli.py clip-map se \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_se \
  --rMATS data/test/clip/ES/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label PTB \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `a3ss` (Alternative 3' Splice Site)
Use for alternative acceptor CLIP enrichment maps.

```bash
python cli.py clip-map a3ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_a3ss \
  --rMATS data/test/clip/A3SS/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label RBFOX2 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `a5ss` (Alternative 5' Splice Site)
Use for alternative donor CLIP enrichment maps.

```bash
python cli.py clip-map a5ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_a5ss \
  --rMATS data/test/clip/A5SS/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label ESRP \
  --sigFDR 0.005 --sigDeltaPSI 0.01 \
  --separate
```

#### `ri` (Retained Intron)
Use for retained-intron CLIP enrichment maps.

```bash
python cli.py clip-map ri \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_ri \
  --rMATS data/test/clip/RI/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label MBNL1 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

#### `mxe` (Mutually Exclusive Exons)
Use for mutually exclusive exon CLIP enrichment maps.

```bash
python cli.py clip-map mxe \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_mxe \
  --rMATS data/test/clip/MXE/test.rMATS.txt \
  --miso NA --up NA --down NA --background NA \
  --label NOVA1 \
  --sigFDR 0.05 --sigDeltaPSI 0.05
```

## Conversion Commands

```bash
python cli.py convert miso --help
```

Examples:

```bash
python cli.py convert miso --event se --in data/test/ESRP.OE.miso_bf --out temp/se.from_miso.rmats.txt
python cli.py convert miso --event a3ss --in path/to/a3ss.miso_bf --out temp/a3ss.from_miso.rmats.txt
python cli.py convert miso --event a5ss --in path/to/a5ss.miso_bf --out temp/a5ss.from_miso.rmats.txt
python cli.py convert miso --event ri --in path/to/ri.miso_bf --out temp/ri.from_miso.rmats.txt
python cli.py convert miso --event mxe --in path/to/mxe.miso_bf --out temp/mxe.from_miso.rmats.txt
```

## Exon Set Commands

```bash
python cli.py exon-sets se --help
```

Example:

```bash
python cli.py exon-sets se \
  --input path/to/common.txt \
  --sample1 SAMPLE1 \
  --sample2 SAMPLE2 \
  --out temp/exon_sets
```

## Output Notes

Motif-map output typically includes:

- `exon/`
- `fasta/`
- `maps/`
- `temp/`
- `log.motifMap.txt`
- `pVal.up.vs.bg.RNAmap.txt`
- `pVal.dn.vs.bg.RNAmap.txt`

CLIP-map output typically includes:

- `exon/`
- `temp/`
- `*.RNAmap.txt`
- `*.pdf` / `*.eps`
- `log.CLIPSeq*.txt`

