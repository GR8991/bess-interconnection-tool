import numpy as np

# ================================
# CONSTANTS
# ================================
MU0 = 4 * np.pi * 1e-7
EPS0 = 8.854e-12


# ================================
# UNIT CONVERSIONS
# ================================
def meter_to_feet(x):
    return x * 3.28084

def km_to_mile(x):
    return x * 0.621371

def ohm_km_to_ohm_mile(r):
    return r / 1.60934


# ================================
# CARSON P–Q APPROX (FULL ANALYTICAL)
# ================================
def carson_PQ(height_ft, rho, freq):
    depth = np.sqrt(rho / (np.pi * freq * MU0)) / 3.28084  # convert meters→feet
    P = np.log(2 * depth / height_ft) - 1
    Q = 1.2533 * (height_ft / depth)
    return P, Q


def carson_self(r_ac_mile, gmr_ft, height_ft, rho, freq):
    P, Q = carson_PQ(height_ft, rho, freq)
    return r_ac_mile + 1j * 0.12134 * (np.log((2*height_ft)/gmr_ft) + P + 1j*Q)


def carson_mutual(D_ft, h1_ft, h2_ft, rho, freq):
    P, Q = carson_PQ(h1_ft + h2_ft, rho, freq)
    De = np.sqrt(D_ft**2 + (h1_ft + h2_ft)**2)
    return 1j * 0.12134 * (np.log(De/D_ft) + P + 1j*Q)


# ================================
# SEQUENCE TRANSFORMATION
# ================================
def seq_transform(Z):
    a = np.exp(1j * 2*np.pi/3)
    T = (1/3)*np.array([
        [1,1,1],
        [1,a,a*a.conjugate()],
        [1,a*a.conjugate(),a]
    ])
    Z012 = T @ Z @ np.linalg.inv(T)
    return Z012[1,1], Z012[0,0]   # return Z1, Z0


# ================================
# MAIN FUNCTION (USED BY your app.py)
# ================================
def compute_overhead_impedance(
    r_ac_km, gmr_m, diameter_m,
    Dab_m, Dbc_m, Dca_m,
    ha_m, hb_m, hc_m,
    rho,
    length,
    unit="mile",
    freq=60
):
    # Unit conversion
    r_ac_mile = ohm_km_to_ohm_mile(r_ac_km)
    gmr_ft = meter_to_feet(gmr_m)

    Dab_ft = meter_to_feet(Dab_m)
    Dbc_ft = meter_to_feet(Dbc_m)
    Dca_ft = meter_to_feet(Dca_m)

    ha_ft = meter_to_feet(ha_m)
    hb_ft = meter_to_feet(hb_m)
    hc_ft = meter_to_feet(hc_m)

    if unit == "mile":
        L_mile = length
    else:
        L_mile = km_to_mile(length)

    # Build Z-phase matrix
    Z = np.zeros((3,3), dtype=complex)

    # Self impedances
    Z[0,0] = carson_self(r_ac_mile, gmr_ft, ha_ft, rho, freq)
    Z[1,1] = carson_self(r_ac_mile, gmr_ft, hb_ft, rho, freq)
    Z[2,2] = carson_self(r_ac_mile, gmr_ft, hc_ft, rho, freq)

    # Mutual
    Z[0,1] = Z[1,0] = carson_mutual(Dab_ft, ha_ft, hb_ft, rho, freq)
    Z[1,2] = Z[2,1] = carson_mutual(Dbc_ft, hb_ft, hc_ft, rho, freq)
    Z[0,2] = Z[2,0] = carson_mutual(Dca_ft, ha_ft, hc_ft, rho, freq)

    # Multiply by line length
    Z = Z * L_mile

    # Sequence extraction
    Z1, Z0 = seq_transform(Z)

    # Shunt (simple model)
    radius_ft = meter_to_feet(diameter_m) / 2
    C = (2 * np.pi * EPS0) / np.log(Dab_ft / radius_ft)
    omega = 2 * np.pi * freq
    B1 = omega * C * L_mile * 5280
    B0 = 0.75 * B1

    # Return EXACT format your app.py expects
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
return Z1_dict, Z0_dict, float(B1), float(B0)


# ================================================================
# PLACEHOLDER FUNCTIONS (fixes import error in app.py)
# ================================================================
def compute_underground_impedance(*args, **kwargs):
    return 0, 0, 0, 0

def compute_multicircuit(*args, **kwargs):
    return {"status": "placeholder - not implemented"}
