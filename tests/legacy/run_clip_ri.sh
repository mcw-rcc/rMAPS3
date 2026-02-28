#!/usr/bin/env bash
#
# Test CLIP-seq RNA map generation for RI events
#

set -e

echo "Testing CLIP-seq RI analysis..."

python cli.py clip-map ri \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_test_ri \
  --rMATS data/test/clip/RI/ri.rMATS.txt \
  --miso NA \
  --up NA \
  --down NA \
  --background NA \
  --label TestRBP

echo "CLIP-seq RI test completed! Check results/clip_test_ri/"


