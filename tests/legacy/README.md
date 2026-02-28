# Legacy Test Scripts

These shell scripts have been **superseded** by the unified Python test scripts:
- `../test_clip.py` - for all CLIP Map tests
- `../test_motif.py` - for all Motif Map tests

## Why They Were Moved

The shell scripts created too much clutter with 14+ individual files for what can be handled by 2 unified Python scripts with command-line options.

## Migration Guide

### Old Shell Scripts to New Python Scripts

**CLIP Map Tests:**
```bash
# OLD: bash tests/run_clip_se.sh
# NEW:
python tests/test_clip.py --event se

# OLD: bash tests/run_all_clip.sh
# NEW:
python tests/test_clip.py
```

**Motif Map Tests:**
```bash
# OLD: bash tests/run_se.sh
# NEW:
python tests/test_motif.py --event se

# OLD: bash tests/run_all_events.sh
# NEW:
python tests/test_motif.py
```

**Multiple Event Types:**
```bash
# NEW: Test multiple events at once
python tests/test_clip.py --event se a5ss a3ss
python tests/test_motif.py --event se a5ss
```

## Compatibility

These legacy scripts are kept for backward compatibility but are **not maintained**. They may be removed in future versions.

For current testing, please use:
- `tests/test_clip.py` 
- `tests/test_motif.py`
- `tests/smoke_cli.py`

See `../README.md` for full documentation.
