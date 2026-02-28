#!/usr/bin/env python3
"""
Unified CLIP Map Test Suite

Usage:
    python tests/test_clip.py            # Test all event types
    python tests/test_clip.py --event se # Test specific event type
    python tests/test_clip.py --list     # List available event types
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime

# CLIP event type configurations
CLIP_EVENTS = {
    "se": {
        "name": "Skipped Exon (SE)",
        "rMATS": "data/test/clip/ES/test.rMATS.txt",
        "peak": "data/test/clip/PIPE-CLIP.Clusters.bed",
        "output": "results/test_clip_se",
    },
    "a5ss": {
        "name": "Alternative 5' Splice Site (A5SS)",
        "rMATS": "data/test/clip/A5SS/test.rMATS.txt",
        "peak": "data/test/clip/PIPE-CLIP.Clusters.bed",
        "output": "results/test_clip_a5ss",
    },
    "a3ss": {
        "name": "Alternative 3' Splice Site (A3SS)",
        "rMATS": "data/test/clip/A3SS/test.rMATS.txt",
        "peak": "data/test/clip/PIPE-CLIP.Clusters.bed",
        "output": "results/test_clip_a3ss",
    },
    "mxe": {
        "name": "Mutually Exclusive Exons (MXE)",
        "rMATS": "data/test/clip/MXE/test.rMATS.txt",
        "peak": "data/test/clip/PIPE-CLIP.Clusters.bed",
        "output": "results/test_clip_mxe",
    },
    "ri": {
        "name": "Retained Intron (RI)",
        "rMATS": "data/test/clip/RI/test.rMATS.txt",
        "peak": "data/test/clip/PIPE-CLIP.Clusters.bed",
        "output": "results/test_clip_ri",
    },
}

# Common parameters for all CLIP tests
COMMON_PARAMS = {
    "label": "TestRBP",
    "miso": "NA",
    "up": "NA",
    "down": "NA",
    "background": "NA",
}

# Match CLI defaults by event type.
EVENT_THRESHOLDS = {
    "se": {"sigFDR": "0.05", "sigDeltaPSI": "0.05"},
    "a3ss": {"sigFDR": "0.05", "sigDeltaPSI": "0.05"},
    "a5ss": {"sigFDR": "0.005", "sigDeltaPSI": "0.01"},
    "mxe": {"sigFDR": "0.05", "sigDeltaPSI": "0.05"},
    "ri": {"sigFDR": "0.05", "sigDeltaPSI": "0.05"},
}


def run_clip_test(event_type, config, verbose=False):
    """Run CLIP test for a specific event type."""
    print(f"\n{'='*70}")
    print(f"Testing CLIP-seq: {config['name']}")
    print(f"{'='*70}")

    thresholds = EVENT_THRESHOLDS[event_type]
    cmd = [
        sys.executable, "cli.py", "clip-map", event_type,
        "--peak", config["peak"],
        "--output", config["output"],
        "--rMATS", config["rMATS"],
        "--miso", COMMON_PARAMS["miso"],
        "--up", COMMON_PARAMS["up"],
        "--down", COMMON_PARAMS["down"],
        "--background", COMMON_PARAMS["background"],
        "--label", COMMON_PARAMS["label"],
        "--sigFDR", thresholds["sigFDR"],
        "--sigDeltaPSI", thresholds["sigDeltaPSI"],
    ]

    if verbose:
        print(f"Command: {' '.join(cmd)}")

    start_time = datetime.now()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent
        )

        elapsed = (datetime.now() - start_time).total_seconds()

        if result.returncode == 0:
            print(f"PASS in {elapsed:.1f}s")
            if verbose and result.stdout:
                print(f"Output: {result.stdout}")
            return True

        print(f"FAIL (exit code {result.returncode})")
        if result.stderr:
            print(f"Error: {result.stderr}")
        if result.stdout:
            print(f"Output: {result.stdout}")
        return False

    except Exception as e:
        print(f"FAIL with exception: {e}")
        return False


def list_events():
    """List all available event types."""
    print("\nAvailable CLIP event types:")
    print("-" * 50)
    for event_type, config in CLIP_EVENTS.items():
        print(f"  {event_type:6s} - {config['name']}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Unified CLIP Map Test Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_clip.py                    # Test all event types
  python tests/test_clip.py --event se         # Test SE only
  python tests/test_clip.py --event se a5ss    # Test SE and A5SS
  python tests/test_clip.py --verbose          # Show detailed output
  python tests/test_clip.py --list             # List event types
        """
    )

    parser.add_argument(
        "--event",
        nargs="+",
        choices=list(CLIP_EVENTS.keys()),
        help="Specific event type(s) to test (default: all)"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available event types and exit"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )

    args = parser.parse_args()

    if args.list:
        list_events()
        return 0

    events_to_test = args.event if args.event else list(CLIP_EVENTS.keys())

    print("\n" + "="*70)
    print("CLIP MAP TEST SUITE")
    print("="*70)
    print(f"Testing {len(events_to_test)} event type(s): {', '.join(events_to_test)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    results = {}
    for event_type in events_to_test:
        config = CLIP_EVENTS[event_type]
        success = run_clip_test(event_type, config, args.verbose)
        results[event_type] = success

    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)

    passed = sum(1 for success in results.values() if success)
    failed = len(results) - passed

    for event_type, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {event_type:6s} - {status}")

    print("-" * 70)
    print(f"Total: {len(results)} tests, {passed} passed, {failed} failed")
    print("="*70)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

