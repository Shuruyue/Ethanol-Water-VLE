"""Tests for output-pack generation."""

from pathlib import Path
import shutil
import unittest

from output_pack import generate_output_pack


class TestOutputPack(unittest.TestCase):
    """Validate practical output export."""

    def test_generate_output_pack(self):
        out_dir = Path(__file__).resolve().parents[1] / "figures" / "_test_pack"
        if out_dir.exists():
            shutil.rmtree(out_dir)

        files = generate_output_pack(
            system_id="ethanol_water_1atm",
            output_dir=out_dir,
            points=30,
            t_ref_celsius=78.2,
        )
        self.assertGreaterEqual(len(files), 8)
        for item in files:
            self.assertTrue(item.exists())

        shutil.rmtree(out_dir)


if __name__ == "__main__":
    unittest.main()
