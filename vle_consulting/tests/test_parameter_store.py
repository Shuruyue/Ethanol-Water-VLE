"""Tests for parameter database loading."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from parameter_store import load_system_parameters


class TestParameterStore(unittest.TestCase):
    """Validate source metadata and NRTL structure."""

    def test_load_default_system(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.assertEqual(system.system_id, "ethanol_water_1atm")
        self.assertGreater(system.pressure_kpa, 100.0)

    def test_nrtl_constants_not_both_zero(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.assertFalse(system.nrtl.a12 == 0.0 and system.nrtl.a21 == 0.0)

    def test_source_metadata_present(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.assertTrue(system.component_1.source.name)
        self.assertTrue(system.component_2.source.citation)
        self.assertTrue(system.nrtl_source.citation)


if __name__ == "__main__":
    unittest.main()
