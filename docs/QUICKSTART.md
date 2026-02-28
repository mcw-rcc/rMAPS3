# Quickstart

## Run CLI Help

```bash
python cli.py --help
python cli.py motif-map --help
```

## Quick SE Run (rMATS Input)

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/example_se \
  --rMATS testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt \
  --miso NA \
  --up NA --down NA --background NA
```

## Quick SE Run (MISO Input)

```bash
python cli.py motif-map se \
  --known-motifs data/knownMotifs.human.mouse.txt \
  --motifs data/ESRP.like.motif.txt \
  --fasta-root genomedata \
  --genome hg19 \
  --output temp/example_se_miso \
  --rMATS NA \
  --miso testData/ESRP.OE.miso_bf \
  --up NA --down NA --background NA
```

## Run Full Event Matrix

```bash
bash tests/run_all_events.sh
```
