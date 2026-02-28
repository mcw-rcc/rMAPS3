#!/usr/bin/env bash
#
# Run all CLIP-seq tests for all event types
#

set -e

echo "=========================================="
echo "Running all CLIP-seq event type tests"
echo "=========================================="

./tests/run_clip_se.sh
./tests/run_clip_a3ss.sh
./tests/run_clip_a5ss.sh
./tests/run_clip_ri.sh
./tests/run_clip_mxe.sh

echo "=========================================="
echo "All CLIP-seq tests completed!"
echo "=========================================="
