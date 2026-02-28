# CLI Usage

## Global Help

```bash
python cli.py --help
python cli.py motif-map --help
python cli.py convert --help
python cli.py exon-sets --help
```

## Motif Map Commands

Event subcommands:

```bash
python cli.py motif-map se --help
python cli.py motif-map a3ss --help
python cli.py motif-map a5ss --help
python cli.py motif-map ri --help
python cli.py motif-map mxe --help
```

Shared required options for motif-map runs:

- `--known-motifs`
- `--fasta-root`
- `--genome`
- `--output`

Common optional options:

- `--motifs` (default: `NA`)
- `--label`
- `--intron`, `--exon`, `--window`, `--step`
- `--sig-fdr`, `--sig-delta-psi`
- `--separate`

### Example: SE from rMATS

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/se_rmats \
  --rMATS testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

### Example: SE from MISO

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/se_miso \
  --rMATS NA \
  --miso testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA
```

### Example: SE from coordinate triplet

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/se_coords \
  --rMATS NA --miso NA \
  --up path/to/up.coord.txt \
  --down path/to/dn.coord.txt \
  --background path/to/bg.coord.txt
```

### Example: A3SS from rMATS

```bash
python cli.py motif-map a3ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/a3ss_rmats \
  --rMATS temp/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

### Example: A5SS from rMATS

```bash
python cli.py motif-map a5ss \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/a5ss_rmats \
  --rMATS temp/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

### Example: RI from rMATS

```bash
python cli.py motif-map ri \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/ri_rmats \
  --rMATS temp/RI.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

### Example: MXE from rMATS

```bash
python cli.py motif-map mxe \
  --known-motifs data/testMotifs.txt \
  --motifs NA \
  --fasta-root genomedata \
  --genome hg19 \
  --output results/mxe_rmats \
  --rMATS temp/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

## Convert Commands

```bash
python cli.py convert miso --help
```

SE conversion:

```bash
python cli.py convert miso --event se --in testData/ESRP.OE.miso_bf --out temp/se.from_miso.rmats.txt
```

A3SS conversion:

```bash
python cli.py convert miso --event a3ss --in path/to/a3ss.miso_bf --out temp/a3ss.from_miso.rmats.txt
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

Motif-map runs typically create:

- `exon/`
- `fasta/`
- `maps/`
- `temp/`
- `log.motifMap.txt`
- `pVal.up.vs.bg.RNAmap.txt`
- `pVal.dn.vs.bg.RNAmap.txt`
