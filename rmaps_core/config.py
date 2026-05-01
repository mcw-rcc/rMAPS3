import os
from pathlib import Path
import warnings

def get_repo_root() -> Path:
    env_repo_root = os.environ.get("RMAPS3_REPO_ROOT")

    if env_repo_root:
        repo_root = Path(env_repo_root).expanduser().resolve()

        if not repo_root.is_dir():
            raise RuntimeError(f"RMAPS3_REPO_ROOT is set but is not a valid directory: {env_repo_root}")

        return repo_root

    warnings.warn(
        "RMAPS3_REPO_ROOT not set; falling back to script location. "
        "This mode is unsupported on HPC systems.",
        RuntimeWarning,
    )

    return Path(__file__).resolve().parents[1]