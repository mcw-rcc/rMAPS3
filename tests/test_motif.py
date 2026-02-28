#!/usr/bin/env python3
"""
Unified Motif Map Test Suite

Usage:
    python tests/test_motif.py            # Test all event types
    python tests/test_motif.py --event se # Test specific event type
    python tests/test_motif.py --list     # List available event types
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime

# Motif event type configurations.
# Uses bundled event-specific rMATS files under data/test/clip/.
MOTIF_EVENTS = {
    "se": {
        "name": "Skipped Exon (SE)",
        "rMATS": "data/test/clip/ES/test.rMATS.txt",
        "output": "results/test_motif_se",
    },
    "a5ss": {
        "name": "Alternative 5' Splice Site (A5SS)",
        "rMATS": "data/test/clip/A5SS/test.rMATS.txt",
        "output": "results/test_motif_a5ss",
    },
    "a3ss": {
        "name": "Alternative 3' Splice Site (A3SS)",
        "rMATS": "data/test/clip/A3SS/test.rMATS.txt",
        "output": "results/test_motif_a3ss",
    },
    "mxe": {
        "name": "Mutually Exclusive Exons (MXE)",
        "rMATS": "data/test/clip/MXE/test.rMATS.txt",
        "output": "results/test_motif_mxe",
    },
    "ri": {
        "name": "Retained Intron (RI)",
        "rMATS": "data/test/clip/RI/test.rMATS.txt",
        "output": "results/test_motif_ri",
    },
}

COMMON_PARAMS = {
    "genome": "hg19",
    "known_motifs": "data/testMotifs.txt",
    "motifs": "NA",
    "miso": "NA",
    "up": "NA",
    "down": "NA",
    "background": "NA",
}


def run_motif_test(event_type, config, fasta_root, verbose=False):
    """Run motif-map test for a specific event type."""
    print(f"\n{'='*70}")
    print(f"Testing Motif Map: {config['name']}")
    print(f"{'='*70}")

    cmd = [
        sys.executable, "cli.py", "motif-map", event_type,
        "--known-motifs", COMMON_PARAMS["known_motifs"],
        "--motifs", COMMON_PARAMS["motifs"],
        "--fasta-root", fasta_root,
        "--genome", COMMON_PARAMS["genome"],
        "--output", config["output"],
        "--rMATS", config["rMATS"],
        "--miso", COMMON_PARAMS["miso"],
        "--up", COMMON_PARAMS["up"],
        "--down", COMMON_PARAMS["down"],
        "--background", COMMON_PARAMS["background"],
    ]

    if verbose:
        print(f"Command: {' '.join(cmd)}")

    start_time = datetime.now()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
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
    print("\nAvailable Motif event types:")
    print("-" * 50)
    for event_type, config in MOTIF_EVENTS.items():
        print(f"  {event_type:6s} - {config['name']}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Unified Motif Map Test Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_motif.py                    # Test all event types
  python tests/test_motif.py --event se         # Test SE only
  python tests/test_motif.py --event se a5ss    # Test SE and A5SS
  python tests/test_motif.py --verbose          # Show detailed output
  python tests/test_motif.py --list             # List event types
        """,
    )

    parser.add_argument(
        "--event",
        nargs="+",
        choices=list(MOTIF_EVENTS.keys()),
        help="Specific event type(s) to test (default: all)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available event types and exit",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output",
    )
    parser.add_argument(
        "--fasta-root",
        default="genomedata",
        help="Path to FASTA root directory (default: genomedata)",
    )
    parser.add_argument(
        "--genome",
        default="hg19",
        help="Genome build to use (default: hg19)",
    )

    args = parser.parse_args()

    if args.list:
        list_events()
        return 0

    COMMON_PARAMS["genome"] = args.genome
    events_to_test = args.event if args.event else list(MOTIF_EVENTS.keys())

    print("\n" + "=" * 70)
    print("MOTIF MAP TEST SUITE")
    print("=" * 70)
    print(f"Testing {len(events_to_test)} event type(s): {', '.join(events_to_test)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"FASTA root: {args.fasta_root}")
    print(f"Genome: {args.genome}")

    results = {}
    for event_type in events_to_test:
        config = MOTIF_EVENTS[event_type]
        success = run_motif_test(event_type, config, args.fasta_root, args.verbose)
        results[event_type] = success

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    passed = sum(1 for success in results.values() if success)
    failed = len(results) - passed

    for event_type, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {event_type:6s} - {status}")

    print("-" * 70)
    print(f"Total: {len(results)} tests, {passed} passed, {failed} failed")
    print("=" * 70)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

