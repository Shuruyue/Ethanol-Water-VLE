"""Shared test configuration for VLE consulting tests."""

import sys
from pathlib import Path

# Add src directory to path so tests can import VLE modules without sys.path hacks
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
