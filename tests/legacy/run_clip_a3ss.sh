#!/usr/bin/env bash
#
# Test CLIP-seq RNA map generation for A3SS events
#

set -e

echo "Testing CLIP-seq A3SS analysis..."

python cli.py clip-map a3ss \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_test_a3ss \
  --rMATS data/test/clip/A3SS/a3ss.rMATS.txt \
  --miso NA \
  --up NA \
  --down NA \
  --background NA \
  --label TestRBP

echo "CLIP-seq A3SS test completed! Check results/clip_test_a3ss/"


