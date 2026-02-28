# CLI Usage

## Global

```bash
python cli.py --help
```

## Motif Map Commands

```bash
python cli.py motif-map --help
python cli.py motif-map se --help
python cli.py motif-map a3ss --help
python cli.py motif-map a5ss --help
python cli.py motif-map ri --help
python cli.py motif-map mxe --help
```

## Convert Commands

```bash
python cli.py convert --help
python cli.py convert miso --help
python cli.py convert miso --event se --in input.miso_bf --out output.rmats.txt
```

## Exon Set Commands

```bash
python cli.py exon-sets --help
python cli.py exon-sets se --help
python cli.py exon-sets se --input path/to/common.txt --sample1 SAMPLE1 --sample2 SAMPLE2 --out temp/exon_sets
```
