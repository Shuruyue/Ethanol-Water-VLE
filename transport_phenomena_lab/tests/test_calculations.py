"""
Unit Tests for Transport Phenomena Lab
單元測試

Run with: python -m unittest tests/test_calculations.py
"""

import unittest
import math
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parents[1]))

from mass_transfer import (
    ETHANOL, WATER,
    antoine_vapor_pressure,
    diffusion_coefficient_fsg,
    mass_transfer_natural,
    mass_transfer_forced,
)

class TestProperties(unittest.TestCase):
    """Test chemical property calculations."""
    
    def test_ethanol_vapor_pressure(self):
        """Verify ethanol vapor pressure at 25°C."""
        # Accepted value ~7.87 kPa (NIST)
        # Calculated using Antoine
        P_sat = antoine_vapor_pressure(25.0, ETHANOL)
        self.assertAlmostEqual(P_sat, 7.83, delta=0.1)
        
    def test_water_vapor_pressure(self):
        """Verify water vapor pressure at 100°C."""
        # Should be ~101.325 kPa
        P_sat = antoine_vapor_pressure(100.0, WATER)
        self.assertAlmostEqual(P_sat, 101.3, delta=0.5)

class TestDiffusion(unittest.TestCase):
    """Test diffusion coefficient calculations."""
    
    def test_fsg_coefficient(self):
        """Verify FSG correlation for Ethanol-Air at 25°C."""
        # Literature value approx 1.18e-5 to 1.35e-5 m2/s
        D_ab = diffusion_coefficient_fsg(25.0, 1.0, ETHANOL)
        self.assertGreater(D_ab, 1.0e-5)
        self.assertLess(D_ab, 1.5e-5)
        
    def test_temperature_dependence(self):
        """Diffusion coefficient should increase with T."""
        D_25 = diffusion_coefficient_fsg(25.0, 1.0, ETHANOL)
        D_50 = diffusion_coefficient_fsg(50.0, 1.0, ETHANOL)
        self.assertTrue(D_50 > D_25)

class TestMassTransfer(unittest.TestCase):
    """Test mass transfer coefficients."""
    
    def test_natural_convection(self):
        """Test natural convection calculation."""
        km, Sh = mass_transfer_natural(25.0, 0.08, ETHANOL)
        self.assertGreater(km, 0)
        self.assertGreater(Sh, 2.0)  # Sh for diffusion only is 2.0
        
    def test_forced_convection(self):
        """Test forced convection calculation."""
        # Forced km should be greater than natural km at 1.5 m/s
        km_nat, _ = mass_transfer_natural(25.0, 0.08, ETHANOL)
        km_for, _ = mass_transfer_forced(25.0, 1.5, 0.08, ETHANOL)
        self.assertGreater(km_for, km_nat)
        
    def test_enhancement_factor(self):
        """Enhancement factor should be reasonable (e.g. > 1.5x at 1.5m/s)."""
        km_nat, _ = mass_transfer_natural(25.0, 0.08, ETHANOL)
        km_for, _ = mass_transfer_forced(25.0, 1.5, 0.08, ETHANOL)
        factor = km_for / km_nat
        self.assertGreater(factor, 1.5)

if __name__ == '__main__':
    unittest.main()
