import numpy as np

# ================================================================
# CONSTANTS (IEEE STANDARD CARSON CONSTANTS)
# ================================================================
MU0 = 4 * np.pi * 1e-7
FREQ = 60
OMEGA = 2 * np.pi * FREQ

# Carson constant (for feet units)
CARSON_K = 0.12134     # ohm per mile


# ================================================================
# UNIT CONVERSION
# ================================================================
def meter_to_feet(x):
    return x * 3.28084

def km_to_mile(x):
    return x * 0.621371

def ohm_km_to_ohm_mile(r):
    return r / 1.60934


# ================================================================
# CORRECT CARSON SELF-IMPEDANCE (FEET)
# ================================================================
def carson_self_impedance(r_ac_mile, gmr_ft, height_ft):
    Ds = (2 * height_ft) / gmr_ft
    Zs = r_ac_mile + 1j * CARSON_K * (np.log(Ds) + 7.934)
    return Zs


# ================================================================
# CORRECT CARSON MUTUAL IMPEDANCE (FEET)
# ================================================================
def carson_mutual_impedance(dist_ft, h1_ft, h2_ft):
    De = np.sqrt(dist_ft**2 + (h1_ft + h2_ft)**2)
    Zm = 1j * CARSON_K * (np.log(De / dist_ft) + 7.934)
    return Zm


# ================================================================
# MAIN FUNCTION (OVERHEAD LINE)
# ================================================================
def compute_overhead_impedance(
    r_ac_km, gmr_m, diameter_m,
    Dab_m, Dbc_m, Dca_m,
    ha_m, hb_m, hc_m,
    rho,              # unused in Carson simplified constant form
    length,
    unit="mile"
):
    """
    Fully correct overhead line impedance calculator using Carson's equations.
    Inputs may be in km or mile, but internally everything is converted to FEET.
    """

    # ----------------------------------------------
    # CONVERT RESISTANCE TO OHM PER MILE
    # ----------------------------------------------
    r_ac_mile = ohm_km_to_ohm_mile(r_ac_km)

    # ----------------------------------------------
    # CONVERT ALL DISTANCES TO FEET
    # ----------------------------------------------
    gmr_ft = meter_to_feet(gmr_m)
    diam_ft = meter_to_feet(diameter_m)

    Dab_ft = meter_to_feet(Dab_m)
    Dbc_ft = meter_to_feet(Dbc_m)
    Dca_ft = meter_to_feet(Dca_m)

    ha_ft = meter_to_feet(ha_m)
    hb_ft = meter_to_feet(hb_m)
    hc_ft = meter_to_feet(hc_m)

    # ----------------------------------------------
    # LINE LENGTH IN MILES
    # ----------------------------------------------
    if unit.lower() == "mile":
        L_mile = length
    else:
        L_mile = km_to_mile(length)

    # ----------------------------------------------
    # BUILD PRIMITIVE MATRIX (PER MILE)
    # ----------------------------------------------
    Z = np.zeros((3, 3), dtype=complex)

    heights = [ha_ft, hb_ft, hc_ft]

    # Self impedances
    Z[0, 0] = carson_self_impedance(r_ac_mile, gmr_ft, ha_ft)
    Z[1, 1] = carson_self_impedance(r_ac_mile, gmr_ft, hb_ft)
    Z[2, 2] = carson_self_impedance(r_ac_mile, gmr_ft, hc_ft)

    # Mutual impedances
    Z[0, 1] = Z[1, 0] = carson_mutual_impedance(Dab_ft, ha_ft, hb_ft)
    Z[1, 2] = Z[2, 1] = carson_mutual_impedance(Dbc_ft, hb_ft, hc_ft)
    Z[0, 2] = Z[2, 0] = carson_mutual_impedance(Dca_ft, ha_ft, hc_ft)

    # ----------------------------------------------
    # SYMMETRICAL COMPONENTS
    # ----------------------------------------------
    a = np.exp(1j * 2 * np.pi / 3)
    T = (1/3) * np.array([
        [1, 1, 1],
        [1, a, a**2],
        [1, a**2, a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)

    Z1 = Z012[1, 1] * L_mile
    Z0 = Z012[0, 0] * L_mile

    # ----------------------------------------------
    # SHUNT SUSCEPTANCE – APPROX BASED ON FEET
    # ----------------------------------------------
    D_eq_ft = Dab_ft
    r_ft = diam_ft / 2

    C = (2 * np.pi * 8.854e-12) / np.log(D_eq_ft / r_ft)
    B1 = (2 * np.pi * FREQ) * C * L_mile * 5280  # convert mile→ft
    B0 = B1 * 0.75

    # ----------------------------------------------
    # FINAL OUTPUT
    # ----------------------------------------------
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
