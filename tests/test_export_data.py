"""Tests for final-build data export."""

from pathlib import Path
import shutil
import unittest

import numpy as np

from export_data import export_full_data_bundle
from pipeline import run_analysis


class TestFinalBuildExport(unittest.TestCase):
    """Validate final export artifacts."""

    def test_export_full_data_bundle(self):
        out_root = Path(__file__).resolve().parents[1] / "final_outputs" / "_test_data"
        if out_root.exists():
            shutil.rmtree(out_root)

        result = run_analysis(points=30, save_plots=False)
        files = export_full_data_bundle(
            result=result,
            t_ref_celsius=78.2,
            points=30,
            pressure_values_kpa=np.linspace(80.0, 120.0, 5),
            output_dir=out_root,
        )

        self.assertGreaterEqual(len(files), 5)
        for item in files:
            self.assertTrue(item.exists())

        shutil.rmtree(out_root)


if __name__ == "__main__":
    unittest.main()
