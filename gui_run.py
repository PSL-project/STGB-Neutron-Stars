"""Launch the desktop GUI from source or from a PyInstaller bundle."""
import sys
from pathlib import Path

if not getattr(sys, "frozen", False):
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from stgb.gui import main

if __name__ == "__main__":
    raise SystemExit(main())
