import numpy as np

# ================================================================
# CONSTANTS (IEEE Standard)
# ================================================================
MU0 = 4 * np.pi * 1e-7
FREQ = 60
OMEGA = 2 * np.pi * FREQ

CARSON_K = 0.12134     # ohm/mile constant used in Carson eq.


# ================================================================
# UNIT CONVERSION
# ================================================================
def meter_to_feet(x):
    return x * 3.28084

def km_to_mile(x):
    return x * 0.621371

def ohm_km_to_ohm_mile(x):
    return x / 1.60934


# ================================================================
# CARSON SELF IMPEDANCE
# ================================================================
def carson_self_impedance(r_ac_mile, gmr_ft, height_ft):
    Ds = (2 * height_ft) / gmr_ft
    return r_ac_mile + 1j * CARSON_K * (np.log(Ds) + 7.934)


# ================================================================
# CARSON MUTUAL IMPEDANCE
# ================================================================
def carson_mutual_impedance(D_ft, h1_ft, h2_ft):
    De = np.sqrt(D_ft**2 + (h1_ft + h2_ft)**2)
    return 1j * CARSON_K * (np.log(De / D_ft) + 7.934)


# ================================================================
# OVERHEAD LINE IMPEDANCE (Carson Full)
# ================================================================
def compute_overhead_impedance(
    r_ac_km, gmr_m, diameter_m,
    Dab_m, Dbc_m, Dca_m,
    ha_m, hb_m, hc_m,
    rho,              # unused but kept for compatibility
    length,
    unit="mile"
):
    # Convert resistance to ohm/mile
    r_ac_mile = ohm_km_to_ohm_mile(r_ac_km)

    # Convert to feet
    gmr_ft = meter_to_feet(gmr_m)
    diam_ft = meter_to_feet(diameter_m)

    Dab_ft = meter_to_feet(Dab_m)
    Dbc_ft = meter_to_feet(Dbc_m)
    Dca_ft = meter_to_feet(Dca_m)

    ha_ft = meter_to_feet(ha_m)
    hb_ft = meter_to_feet(hb_m)
    hc_ft = meter_to_feet(hc_m)

    # Convert line length
    if unit.lower() == "mile":
        L_mile = length
    else:
        L_mile = km_to_mile(length)

    # Primitive matrix per mile
    Z = np.zeros((3,3), dtype=complex)

    # Self
    Z[0,0] = carson_self_impedance(r_ac_mile, gmr_ft, ha_ft)
    Z[1,1] = carson_self_impedance(r_ac_mile, gmr_ft, hb_ft)
    Z[2,2] = carson_self_impedance(r_ac_mile, gmr_ft, hc_ft)

    # Mutual
    Z[0,1] = Z[1,0] = carson_mutual_impedance(Dab_ft, ha_ft, hb_ft)
    Z[1,2] = Z[2,1] = carson_mutual_impedance(Dbc_ft, hb_ft, hc_ft)
    Z[0,2] = Z[2,0] = carson_mutual_impedance(Dca_ft, ha_ft, hc_ft)

    # Symmetrical components
    a = np.exp(1j * 2*np.pi/3)
    T = (1/3)*np.array([
        [1,1,1],
        [1,a,a*a.conjugate()],
        [1,a*a.conjugate(),a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)

    Z1 = Z012[1,1] * L_mile
    Z0 = Z012[0,0] * L_mile

    # Shunt susceptance (simple approximation)
    D_eq_ft = Dab_ft
    r_ft = diam_ft / 2
    C = (2 * np.pi * 8.854e-12) / np.log(D_eq_ft / r_ft)
    B1 = (2 * np.pi * FREQ) * C * L_mile * 5280
    B0 = B1 * 0.75

    return (
        {
            "R1_ohm": float(np.real(Z1)),
            "X1_ohm": float(np.imag(Z1)),
            "Z1_mag": float(np.abs(Z1)),
            "Z1_angle_deg": float(np.angle(Z1, deg=True))
        },
        {
            "R0_ohm": float(np.real(Z0)),
            "X0_ohm": float(np.imag(Z0)),
            "Z0_mag": float(np.abs(Z0)),
            "Z0_angle_deg": float(np.angle(Z0, deg=True))
        },
        float(B1),
        float(B0)
    )


# ================================================================
# UNDERGROUND CABLE IMPEDANCE (Basic)
# ================================================================
def compute_underground_impedance(
    r_ac_km, gmr_m, diameter_m,
    spacing_m, bonding, length_km
):
    # Simple underground model (not Carson)
    r_per_m = r_ac_km / 1000
    Zs = r_per_m + 1j * OMEGA * MU0 / (2*np.pi) * np.log(2*spacing_m / gmr_m)
    Zm = 1j * OMEGA * MU0 / (2*np.pi) * np.log(spacing_m / diameter_m)

    Z = np.array([
        [Zs, Zm, Zm],
        [Zm, Zs, Zm],
        [Zm, Zm, Zs]
    ])

    a = np.exp(1j * 2*np.pi/3)
    T = (1/3)*np.array([
        [1,1,1],
        [1,a,a*a.conjugate()],
        [1,a*a.conjugate(),a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)

    Z1 = Z012[1,1] * (length_km*1000)
    Z0 = Z012[0,0] * (length_km*1000)

    return Z1, Z0, 0.0, 0.0  # B not computed here


# ================================================================
# MULTI-CIRCUIT (Placeholder)
# ================================================================
def compute_multicircuit(n_circuits, include_ohgw, frequency):
    return {
        "Z1": (0.1/n_circuits) * np.exp(1j*75*np.pi/180),
        "Z0": (0.3/n_circuits) * np.exp(1j*85*np.pi/180),
        "GroundWireIncluded": include_ohgw
    }
