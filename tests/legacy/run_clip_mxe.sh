#!/usr/bin/env bash
#
# Test CLIP-seq RNA map generation for MXE events
#

set -e

echo "Testing CLIP-seq MXE analysis..."

python cli.py clip-map mxe \
  --peak data/test/clip/PIPE-CLIP.Clusters.bed \
  --output results/clip_test_mxe \
  --rMATS data/test/clip/MXE/mxe.rMATS.txt \
  --miso NA \
  --up NA \
  --down NA \
  --background NA \
  --label TestRBP

echo "CLIP-seq MXE test completed! Check results/clip_test_mxe/"


