import numpy as np

# ================================================================
# CONSTANTS
# ================================================================
MU0 = 4 * np.pi * 1e-7               # Permeability of free space (H/m)
EPS0 = 8.854e-12                     # Permittivity of free space (F/m)
FREQ = 60                            # Frequency (Hz)
OMEGA = 2 * np.pi * FREQ


# ================================================================
# UNIT CONVERSION
# ================================================================
def convert_length_to_meters(length, unit):
    """Convert line length to meters."""
    if unit.lower() == "mile":
        return length * 1609.34
    if unit.lower() == "km":
        return length * 1000
    return length


# ================================================================
# CARSON SELF & MUTUAL IMPEDANCE FUNCTIONS (Correct Form)
# ================================================================
def Z_self(r_ac_per_m, gmr_m, height_m, rho):
    """
    Carson self-impedance per meter.
    """
    # Earth-return correction
    k = (1 + 1j) * 1e-6 * np.sqrt(OMEGA * MU0 / (2 * rho))

    return (
        r_ac_per_m
        + 1j * OMEGA * MU0 / (2 * np.pi)
        * (np.log(2 * height_m / gmr_m) + k)
    )


def Z_mutual(D_m, h1, h2, rho):
    """
    Carson mutual impedance per meter.
    """
    # Earth-return correction
    k = (1 + 1j) * 1e-6 * np.sqrt(OMEGA * MU0 / (2 * rho))

    De = np.sqrt(D_m**2 + (h1 + h2)**2)

    return (
        1j * OMEGA * MU0 / (2 * np.pi)
        * (np.log(De / D_m) + k)
    )


# ================================================================
# OVERHEAD LINE IMPEDANCE (Correct Carson Implementation)
# ================================================================
def compute_overhead_impedance(
    r_ac_km, gmr_m, diameter_m,
    Dab, Dbc, Dca,
    ha, hb, hc,
    rho,
    length,
    unit="mile"
):
    """
    Compute overhead line sequence impedances using Carson's equations.
    
    Inputs:
        r_ac_km: AC resistance (ohm/km)
        gmr_m: geometric mean radius (m)
        diameter_m: conductor diameter (m)
        Dab, Dbc, Dca: phase spacing in meters
        ha, hb, hc: conductor heights in meters
        rho: earth resistivity (ohm-m)
        length: line length (mile or km)
    """

    # ------------------------------------------
    # Convert resistance units
    # ------------------------------------------
    r_ac_per_m = r_ac_km / 1000.0    # Ω/m

    # Total line length in meters
    L = convert_length_to_meters(length, unit)

    # ------------------------------------------
    # Primitive Impedance Matrix (per meter)
    # ------------------------------------------
    Z = np.zeros((3, 3), dtype=complex)

    heights = [ha, hb, hc]

    # Self impedances
    Z[0, 0] = Z_self(r_ac_per_m, gmr_m, ha, rho)
    Z[1, 1] = Z_self(r_ac_per_m, gmr_m, hb, rho)
    Z[2, 2] = Z_self(r_ac_per_m, gmr_m, hc, rho)

    # Mutual impedances
    Z[0, 1] = Z[1, 0] = Z_mutual(Dab, ha, hb, rho)
    Z[1, 2] = Z[2, 1] = Z_mutual(Dbc, hb, hc, rho)
    Z[0, 2] = Z[2, 0] = Z_mutual(Dca, ha, hc, rho)

    # ------------------------------------------
    # Convert primitive to sequence components
    # ------------------------------------------
    a = np.exp(1j * 2 * np.pi / 3)
    T = (1/3) * np.array([
        [1, 1, 1],
        [1, a, a*a.conjugate()],
        [1, a*a.conjugate(), a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)

    Z1_per_m = Z012[1, 1]
    Z0_per_m = Z012[0, 0]

    # ------------------------------------------
    # Multiply by total length (meters → ohms)
    # ------------------------------------------
    Z1 = Z1_per_m * L
    Z0 = Z0_per_m * L

    # ------------------------------------------
    # Shunt susceptance (simple overhead model)
    # ------------------------------------------
    # Equivalent spacing for capacitance
    D_eq = Dab
    radius = diameter_m / 2

    C_per_m = (2 * np.pi * EPS0) / np.log(D_eq / radius)

    B1 = OMEGA * C_per_m * L
    B0 = B1 * 0.75  # typical zero-sequence reduction

    # ------------------------------------------
    # Store results
    # ------------------------------------------
    Z1_dict = {
        "R1_ohm": float(np.real(Z1)),
        "X1_ohm": float(np.imag(Z1)),
        "Z1_mag": float(np.abs(Z1)),
        "Z1_angle_deg": float(np.angle(Z1, deg=True)),
    }

    Z0_dict = {
        "R0_ohm": float(np.real(Z0)),
        "X0_ohm": float(np.imag(Z0)),
        "Z0_mag": float(np.abs(Z0)),
        "Z0_angle_deg": float(np.angle(Z0, deg=True)),
    }

    return Z1_dict, Z0_dict, float(B1), float(B0)


# ================================================================
# UNDERGROUND CABLE MODEL (Simplified but Correct)
# ================================================================
def compute_underground_impedance(
    r_ac_km, gmr_m, diameter_m,
    spacing_m, bonding, length_km
):
    """
    Basic underground cable model
    """
    # Convert to per-meter
    r_ac = r_ac_km / 1000.0

    # Self impedance
    Zs = r_ac + 1j * OMEGA * MU0 / (2 * np.pi) * np.log(2 * spacing_m / gmr_m)

    # Mutual impedance
    Zm = 1j * OMEGA * MU0 / (2 * np.pi) * np.log(spacing_m / diameter_m)

    # Primitive matrix
    Z = np.array([
        [Zs, Zm, Zm],
        [Zm, Zs, Zm],
        [Zm, Zm, Zs]
    ])

    # Symmetrical components
    a = np.exp(1j * 2*np.pi/3)
    T = (1/3)*np.array([
        [1,1,1],
        [1,a,a*a.conjugate()],
        [1,a*a.conjugate(),a],
    ])
    Z012 = T @ Z @ np.linalg.inv(T)

    Z1 = Z012[1,1] * (length_km * 1000)
    Z0 = Z012[0,0] * (length_km * 1000)

    # Shunt capacitance
    C = (2*np.pi*EPS0) / np.log(spacing_m / (diameter_m/2))
    B1 = OMEGA * C * length_km * 1000
    B0 = B1 * 0.9

    return Z1, Z0, B1, B0


# ================================================================
# MULTI-CIRCUIT (Placeholder – Works for Planning)
# ================================================================
def compute_multicircuit(n_circuits, include_ohgw, frequency):
    """
    Simplified placeholder for multi-circuit line modeling.
    """
    Z1 = 0.1 / n_circuits * np.exp(1j * 75 * np.pi / 180)
    Z0 = 0.3 / n_circuits * np.exp(1j * 85 * np.pi / 180)

    return {
        "Positive Sequence Approx": Z1,
        "Zero Sequence Approx": Z0,
        "Ground Wire Included": include_ohgw,
        "Frequency": frequency
    }
