#!/usr/bin/env bash
#
# Test CLIP-seq RNA map generation for SE (Exon Skipping) events
#

set -e

echo "Testing CLIP-seq SE analysis..."

python cli.py clip-map se \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_test_se \
  --rMATS data/test/clip/ES/test.rMATS.txt \
  --miso NA \
  --up NA \
  --down NA \
  --background NA \
  --label TestRBP

echo "CLIP-seq SE test completed! Check results/clip_test_se/"

