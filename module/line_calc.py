import numpy as np

# ================================================================
# CONSTANTS
# ================================================================
MU0 = 4 * np.pi * 1e-7               # Permeability of free space (H/m)
EPS0 = 8.854e-12                     # Permittivity of free space (F/m)
FREQ = 60                            # Default frequency
OMEGA = 2 * np.pi * FREQ


# ================================================================
# HELPER FUNCTIONS
# ================================================================
def convert_length(length, unit):
    """Convert line length to meters."""
    if unit.lower() == "km":
        return length * 1000
    elif unit.lower() == "mile":
        return length * 1609.34
    else:
        return length


def carson_self_impedance(r_ac, gmr, height, rho):
    """
    Self-impedance using Carson's approximation.
    r_ac: AC resistance per km
    gmr: geometric mean radius (m)
    height: conductor height (m)
    rho: earth resistivity (Ω·m)
    """
    earth_term = (1 + 1j) * 0.000001 * np.sqrt(OMEGA * MU0 / (2 * rho))
    Zs = (r_ac + 
          1j * OMEGA * MU0 / (2 * np.pi) * (np.log(2 * height / gmr) + earth_term))
    return Zs


def carson_mutual_impedance(D, h1, h2, rho):
    """
    Mutual impedance between two overhead conductors.
    """
    De = np.sqrt(D**2 + (h1 + h2)**2)
    earth_term = (1 + 1j) * 0.000001 * np.sqrt(OMEGA * MU0 / (2 * rho))
    Zm = 1j * OMEGA * MU0 / (2 * np.pi) * (np.log(De / D) + earth_term)
    return Zm


def sequence_impedance(Z):
    """Convert primitive Z-matrix → sequence impedances (Z0, Z1)."""
    a = np.exp(1j * 2*np.pi/3)
    T = (1/3) * np.array([
        [1,    1,    1],
        [1,    a,    a**2],
        [1,    a**2, a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)
    Z0 = Z012[0, 0]
    Z1 = Z012[1, 1]
    return Z1, Z0


def shunt_capacitance(diameter, spacing, eps_r=1):
    """
    Approximate line capacitance per phase for overhead lines.
    Returns susceptance B (siemens per meter).
    """
    D_eq = spacing
    r = diameter / 2
    C = (2 * np.pi * EPS0 * eps_r) / np.log(D_eq / r)
    B = OMEGA * C
    return B


# ================================================================
# OVERHEAD LINE IMPEDANCE (CARSON)
# ================================================================
def compute_overhead_impedance(
    r_ac, gmr, diameter,
    Dab, Dbc, Dca,
    ha, hb, hc,
    rho,
    length,
    unit="km"
):
    """
    Returns:
        Z1_dict = {"R1_ohm": ..., "X1_ohm": ..., "Z1_pu": ...}
        Z0_dict = {"R0_ohm": ..., "X0_ohm": ..., "Z0_pu": ...}
        B1, B0 (shunt susceptance)
    """

    L = convert_length(length, unit)

    # Primitive Z matrix
    Z = np.zeros((3, 3), dtype=complex)

    heights = [ha, hb, hc]

    # Self impedances
    for i, h in enumerate(heights):
        Z[i, i] = carson_self_impedance(r_ac, gmr, h, rho)

    # Mutual impedances
    Z[0, 1] = Z[1, 0] = carson_mutual_impedance(Dab, ha, hb, rho)
    Z[1, 2] = Z[2, 1] = carson_mutual_impedance(Dbc, hb, hc, rho)
    Z[0, 2] = Z[2, 0] = carson_mutual_impedance(Dca, ha, hc, rho)

    # Convert to sequence impedances
    Z1, Z0 = sequence_impedance(Z)

    # Convert per-meter to full-line impedance
    Z1_total = Z1 * L
    Z0_total = Z0 * L

    # Shunt susceptance approximation
    B1 = shunt_capacitance(diameter, Dab) * L
    B0 = B1  # Simplified model for overhead lines

    Z1_dict = {
        "R1_ohm": np.real(Z1_total),
        "X1_ohm": np.imag(Z1_total),
        "Z1_mag": np.abs(Z1_total),
        "Z1_angle_deg": np.angle(Z1_total, deg=True)
    }

    Z0_dict = {
        "R0_ohm": np.real(Z0_total),
        "X0_ohm": np.imag(Z0_total),
        "Z0_mag": np.abs(Z0_total),
        "Z0_angle_deg": np.angle(Z0_total, deg=True)
    }

    return Z1_dict, Z0_dict, B1, B0


# ================================================================
# UNDERGROUND CABLE IMPEDANCE
# ================================================================
def compute_underground_impedance(
    r_ac, gmr, diameter, spacing, bonding, length_km
):
    """
    Simplified underground cable impedance model.
    """

    omega = OMEGA

    # Self impedance
    Zs = r_ac + 1j * omega * MU0 / (2 * np.pi) * np.log(2 * spacing / gmr)

    # Mutual impedance (approx)
    Zm = 1j * omega * MU0 / (2 * np.pi) * np.log(spacing / diameter)

    # Build primitive Z
    Z = np.array([
        [Zs, Zm, Zm],
        [Zm, Zs, Zm],
        [Zm, Zm, Zs]
    ])

    Z1, Z0 = sequence_impedance(Z)

    # Total impedance
    L = length_km * 1000
    Z1_total = Z1 * L
    Z0_total = Z0 * L

    # Approximate capacitance
    B1 = shunt_capacitance(diameter, spacing) * L
    B0 = B1

    return Z1_total, Z0_total, B1, B0


# ================================================================
# MULTI-CIRCUIT OVERHEAD LINE MODEL
# ================================================================
def compute_multicircuit(n_circuits, include_ohgw, frequency):
    """
    Simplified multi-circuit line impedance model.
    Not a full Carson solution but suitable for planning-level studies.
    """

    base_impedance = 0.05 * n_circuits  # placeholder

    result = {
        "Positive Sequence Approx (Ω)": base_impedance * np.exp(1j * 70 * np.pi / 180),
        "Zero Sequence Approx (Ω)": base_impedance * 3 * np.exp(1j * 85 * np.pi / 180),
        "Overhead Ground Wire Included": include_ohgw,
        "Frequency (Hz)": frequency
    }

    return result

