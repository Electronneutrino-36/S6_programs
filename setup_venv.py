#!/usr/bin/env python3
"""Create a local virtual environment and install the required packages."""
import subprocess
import sys
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV_DIR = ROOT / ".venv"
REQUIREMENTS = ROOT / "requirements.txt"


def venv_python() -> Path:
    if sys.platform == "win32":
        return VENV_DIR / "Scripts" / "python.exe"
    return VENV_DIR / "bin" / "python"


def main() -> None:
    if not REQUIREMENTS.exists():
        sys.exit(f"Missing {REQUIREMENTS.name}")

    if venv_python().exists():
        print("Virtual environment already exists, reusing it")
    else:
        print("Creating virtual environment...")
        venv.create(VENV_DIR, with_pip=True)

    print("Installing necessary packages...")
    subprocess.check_call(
        [str(venv_python()), "-m", "pip", "install", "-r", str(REQUIREMENTS)]
    )

    activate = r".venv\Scripts\activate" if sys.platform == "win32" else "source .venv/bin/activate"
    print("\nFinished. To use the environment, run:")
    print(f"  {activate}")


if __name__ == "__main__":
    main()