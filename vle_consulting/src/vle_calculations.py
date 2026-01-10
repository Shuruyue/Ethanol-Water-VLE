import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from thermo.interaction_parameters import IPDB
from thermo.nrtl import NRTL
from chemicals import vapor_pressure

# ============================================================
# Ethanol(1) – Water(2) Thermodynamics at 1 atm (101.325 kPa)
# Models: Raoult (Ideal) / NRTL (Non-ideal)
# ============================================================

# Constants
R = 8.314462618  # J/mol/K
P_SYS = 101325.0  # System pressure in Pa (1 atm)
CAS_ETHANOL = '64-17-5'
CAS_WATER = '7732-18-5'

def get_vapor_pressure(cas, T):
    """Calculate vapor pressure (Pa) at Temperature T (K) using Antoine equation from chemicals library."""
    # We use the default method which usually picks the best available correlation
    # Antoine coefficients are handled internally by chemicals/thermo
    try:
        # P_vap is returned in Pa
        return vapor_pressure.vapor_pressure(CASRN=cas, T=T)
    except Exception as e:
        print(f"Error calculating vapor pressure for {cas} at {T}K: {e}")
        return None

def get_nrtl_model(T, xs):
    """
    Create an NRTL model instance for Ethanol-Water at given T and composition xs.
    Uses ChemSep parameters from IPDB.
    """
    # Get parameters from IPDB (Interaction Parameter Database)
    # Note: IPDB returns parameters in the form suitable for the NRTL class
    # For ChemSep NRTL, we need alpha (alphaij) and tau parameters. 
    # ChemSep provides 'bij' parameters which relate to tau via tau_ij = bij / T (if aij=0)
    # However, thermo's NRTL class can take tau_coeffs directly if we format them correctly.
    
    # Actually, let's use the helper method from the test file which is cleaner
    # IPDB.get_ip_asymmetric_matrix returns the matrix for a specific property
    
    # alpha_ij (non-randomness parameter)
    alphas = IPDB.get_ip_asymmetric_matrix('ChemSep NRTL', [CAS_ETHANOL, CAS_WATER], 'alphaij')
    
    # b_ij (energy parameter in K, such that tau_ij = b_ij/T)
    # The database key is 'bij' for ChemSep NRTL
    bijs = IPDB.get_ip_asymmetric_matrix('ChemSep NRTL', [CAS_ETHANOL, CAS_WATER], 'bij')
    
    # We need to construct tau_coeffs in the format expected by NRTL class.
    # The NRTL class expects tau_coeffs or explicit tau_as, tau_bs, etc.
    # tau_ij = A_ij + B_ij/T + ...
    # So we can pass bijs as tau_bs.
    
    # Construct NRTL object
    # We pass bijs as tau_bs to let the model handle temperature dependence: tau = b/T
    model = NRTL(T=T, xs=xs, alpha_cs=alphas, tau_bs=bijs)
    return model

def solve_bubble_point(x1, model_type='nrtl'):
    """
    Solve for bubble point temperature and vapor composition.
    x1: mole fraction of ethanol in liquid
    model_type: 'ideal' or 'nrtl'
    """
    x2 = 1.0 - x1
    xs = [x1, x2]
    
    def objective(T):
        # Calculate P_sat for both components
        Psat1 = get_vapor_pressure(CAS_ETHANOL, T)
        Psat2 = get_vapor_pressure(CAS_WATER, T)
        
        # Guard against vapor pressure calculation failures
        # Return a penalty value to guide brentq away from invalid temperatures
        if Psat1 is None or Psat2 is None:
            # Return signed residual to help brentq bracketing:
            # High T (>350K) → positive residual → search lower
            # Low T (<350K) → negative residual → search higher
            return 1e9 if T > 350 else -1e9
        
        if model_type == 'ideal':
            gamma1, gamma2 = 1.0, 1.0
        else:
            # NRTL model
            model = get_nrtl_model(T, xs)
            gamma1, gamma2 = model.gammas()
            
        # Bubble point equation: P_sys = x1*gamma1*Psat1 + x2*gamma2*Psat2
        P_calc = x1 * gamma1 * Psat1 + x2 * gamma2 * Psat2
        return P_calc - P_SYS

    # Solve for T (expecting boiling points between ~78C and ~100C)
    try:
        T_boil = brentq(objective, 340, 380) # 67C to 107C range
    except ValueError:
        try:
            T_boil = brentq(objective, 300, 400) # Wider range just in case
        except ValueError:
             print(f"Could not solve T_boil for x1={x1}")
             return None, None

    # Calculate vapor phase composition y
    Psat1 = get_vapor_pressure(CAS_ETHANOL, T_boil)
    Psat2 = get_vapor_pressure(CAS_WATER, T_boil)
    
    if model_type == 'ideal':
        gamma1, gamma2 = 1.0, 1.0
    else:
        model = get_nrtl_model(T_boil, xs)
        gamma1, gamma2 = model.gammas()
        
    y1 = (x1 * gamma1 * Psat1) / P_SYS
    return T_boil, y1

def calculate_excess_properties(x1_vals):
    """
    Calculate G^E and H^E (approx approx Delta_H_mix) for NRTL model along the bubble point curve.
    """
    GE_vals = []
    HE_vals = []
    
    for x1 in x1_vals:
        T, y1 = solve_bubble_point(x1, model_type='nrtl')
        xs = [x1, 1.0-x1]
        model = get_nrtl_model(T, xs)
        
        # Excess Gibbs Energy
        # GE = RT * sum(xi * ln(gamma_i))
        # The NRTL class has a GE() method but it returns specific GE (J/mol) ? Let's check docs or source.
        # Based on source code reading: model.GE() returns excess gibbs energy.
        ge = model.GE()
        GE_vals.append(ge)
        
        # Excess Enthalpy
        # H^E = -T^2 * d(G^E/T)/dT
        # The NRTL class has a dGE_dT() method.
        # G^E = H^E - T*S^E  => H^E = G^E + T*S^E
        # S^E = -dG^E/dT
        # So H^E = G^E - T*(dG^E/dT)
        dge_dt = model.dGE_dT()
        he = ge - T * dge_dt
        HE_vals.append(he)
        
    return np.array(GE_vals), np.array(HE_vals)

def main():
    print("Starting VLE calculations for Ethanol-Water system at 1 atm...")
    
    # Composition grid
    x_vals = np.linspace(0.001, 0.999, 50)
    
    # Storage for results
    T_ideal = []
    y_ideal = []
    T_nrtl = []
    y_nrtl = []
    
    print("Calculating Ideal and NRTL profiles...")
    for x in x_vals:
        # Ideal
        Ti, yi = solve_bubble_point(x, model_type='ideal')
        T_ideal.append(Ti - 273.15) # Convert to C
        y_ideal.append(yi)
        
        # NRTL
        Tn, yn = solve_bubble_point(x, model_type='nrtl')
        T_nrtl.append(Tn - 273.15) # Convert to C
        y_nrtl.append(yn)
        
    # Calculate excess properties for NRTL
    print("Calculating excess properties...")
    GE_vals, HE_vals = calculate_excess_properties(x_vals)
    
    # ==========================
    # Plotting
    # ==========================
    
    # 1. T-x-y Diagram
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, T_ideal, 'b--', label='Ideal Liquid (Tx)')
    plt.plot(y_ideal, T_ideal, 'b:', label='Ideal Vapor (Ty)')
    plt.plot(x_vals, T_nrtl, 'r-', linewidth=2, label='NRTL Liquid (Tx)')
    plt.plot(y_nrtl, T_nrtl, 'r-.', linewidth=2, label='NRTL Vapor (Ty)')
    plt.xlabel('Mole Fraction Ethanol')
    plt.ylabel('Temperature (°C)')
    plt.title('T-x-y Diagram: Ethanol-Water at 1 atm')
    plt.legend()
    plt.grid(True)
    plt.savefig('TxY_Ethanol_Water_1atm.png', dpi=300)
    print("Saved TxY_Ethanol_Water_1atm.png")
    
    # 2. y-x Diagram
    plt.figure(figsize=(8, 8))
    plt.plot([0, 1], [0, 1], 'k-', alpha=0.3, label='Reference (y=x)')
    plt.plot(x_vals, y_ideal, 'b--', label='Ideal (Raoult)')
    plt.plot(x_vals, y_nrtl, 'r-', linewidth=2, label='NRTL')
    plt.xlabel('Liquid Mole Fraction Ethanol (x₁)')
    plt.ylabel('Vapor Mole Fraction Ethanol (y₁)')
    plt.title('y-x Diagram: Ethanol-Water at 1 atm')
    plt.legend()
    plt.grid(True)
    plt.savefig('Yx_Ethanol_Water_1atm.png', dpi=300)
    print("Saved Yx_Ethanol_Water_1atm.png")
    
    # 3. Excess Gibbs Energy
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, GE_vals, 'k-', linewidth=2)
    plt.xlabel('Mole Fraction Ethanol (x₁)')
    plt.ylabel('Excess Gibbs Energy Gᴱ (J/mol)')
    plt.title('NRTL Excess Gibbs Energy vs Composition')
    plt.grid(True)
    plt.savefig('GE_vs_x_Ethanol_Water_1atm.png', dpi=300)
    print("Saved GE_vs_x_Ethanol_Water_1atm.png")
    
    # 4. Excess Enthalpy (Heat of Mixing)
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, HE_vals, 'm-', linewidth=2)
    plt.xlabel('Mole Fraction Ethanol (x₁)')
    plt.ylabel('Excess Enthalpy Hᴱ ≈ ΔH_mix (J/mol)')
    plt.title('NRTL Excess Enthalpy (Heat of Mixing) vs Composition')
    plt.grid(True)
    # Add zero line for reference
    plt.axhline(0, color='k', linestyle='-', alpha=0.3)
    plt.savefig('Hmix_vs_x_Ethanol_Water_1atm.png', dpi=300)
    print("Saved Hmix_vs_x_Ethanol_Water_1atm.png")
    
    print("All calculations and plots completed.")

if __name__ == "__main__":
    main()
