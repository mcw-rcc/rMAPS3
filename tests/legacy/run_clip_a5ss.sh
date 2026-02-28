#!/usr/bin/env bash
#
# Test CLIP-seq RNA map generation for A5SS events
#

set -e

echo "Testing CLIP-seq A5SS analysis..."

python cli.py clip-map a5ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_test_a5ss \
  --rMATS data/test/clip/A5SS/test.rMATS.txt \
  --miso NA \
  --up NA \
  --down NA \
  --background NA \
  --label TestRBP \
  --sigFDR 0.005 \
  --sigDeltaPSI 0.01 \
  --separate

echo "CLIP-seq A5SS test completed! Check results/clip_test_a5ss/"

