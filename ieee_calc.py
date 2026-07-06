"""
IEEE/ANSI interconnection impedance engine.

Two calculators:
  * Section 3 - overhead tie line   -> modified Carson's equations
  * Section 5 - MV collector system -> WECC single-machine aggregation

All functions are pure Python (no Streamlit) so they can be unit-tested.
Distances internally in feet, impedances in ohm/mile (line) or ohm (collector).
"""
import math

EPS0 = 8.854e-12          # F/m
MU0 = 4 * math.pi * 1e-7  # H/m
MI = 1609.34              # m per mile
FT = 0.3048               # m per ft

# --------------------------------------------------------------------------
# CONDUCTOR LIBRARY  (overhead, ACSR)
#   OD_in  : outside diameter (in)
#   gmr_ft : geometric mean radius (ft)
#   rac    : 60 Hz AC resistance (ohm/mile) at 75 C  -- OVERRIDE from datasheet
# Geometry (OD, GMR) is stable; resistance varies with temperature/standard.
# --------------------------------------------------------------------------
CONDUCTORS = {
    "ACSR 266.8 Partridge (26/7)": {"OD_in": 0.642, "gmr_ft": 0.0217, "rac": 0.3792},
    "ACSR 336.4 Linnet (26/7)":    {"OD_in": 0.721, "gmr_ft": 0.0243, "rac": 0.3070},
    "ACSR 397.5 Ibis (26/7)":      {"OD_in": 0.783, "gmr_ft": 0.0265, "rac": 0.2590},
    "ACSR 477 Hawk (26/7)":        {"OD_in": 0.858, "gmr_ft": 0.0289, "rac": 0.2160},
    "ACSR 556.5 Dove (26/7)":      {"OD_in": 0.927, "gmr_ft": 0.0313, "rac": 0.1859},
    "ACSR 636 Grosbeak (26/7)":    {"OD_in": 0.990, "gmr_ft": 0.0335, "rac": 0.1618},
    "ACSR 795 Drake (26/7)":       {"OD_in": 1.108, "gmr_ft": 0.0373, "rac": 0.1284},
    "ACSR 954 Rail (45/7)":        {"OD_in": 1.165, "gmr_ft": 0.0386, "rac": 0.1128},
    "ACSR 954 Cardinal (54/7)":    {"OD_in": 1.196, "gmr_ft": 0.0402, "rac": 0.1099},
    "ACSR 1272 Bittern (45/7)":    {"OD_in": 1.345, "gmr_ft": 0.0445, "rac": 0.0849},
    "ACSR 1590 Falcon (54/19)":    {"OD_in": 1.545, "gmr_ft": 0.0523, "rac": 0.0684},
}

# --------------------------------------------------------------------------
# MV CABLE LIBRARY  (single-core shielded XLPE, per phase)
#   R1,X1,R0,X0 : ohm/km   ;  C1 : uF/km   (C0 = C1 for shielded cable)
# Representative planning values -- OVERRIDE with IEC 60287 / manufacturer data.
# --------------------------------------------------------------------------
CABLES = {
    "Al XLPE 1/0 AWG 35 kV":   {"R1": 0.529, "X1": 0.128, "R0": 0.640, "X0": 0.170, "C1": 0.14},
    "Al XLPE 4/0 AWG 35 kV":   {"R1": 0.262, "X1": 0.116, "R0": 0.360, "X0": 0.150, "C1": 0.18},
    "Al XLPE 250 kcmil 35 kV": {"R1": 0.224, "X1": 0.111, "R0": 0.315, "X0": 0.145, "C1": 0.19},
    "Al XLPE 350 kcmil 35 kV": {"R1": 0.161, "X1": 0.106, "R0": 0.250, "X0": 0.140, "C1": 0.22},
    "Al XLPE 500 kcmil 35 kV": {"R1": 0.113, "X1": 0.101, "R0": 0.195, "X0": 0.135, "C1": 0.24},
    "Al XLPE 750 kcmil 35 kV": {"R1": 0.095, "X1": 0.095, "R0": 0.165, "X0": 0.130, "C1": 0.27},
    "Al XLPE 1000 kcmil 35 kV":{"R1": 0.075, "X1": 0.091, "R0": 0.140, "X0": 0.126, "C1": 0.30},
    "Cu XLPE 500 kcmil 35 kV": {"R1": 0.069, "X1": 0.099, "R0": 0.130, "X0": 0.133, "C1": 0.25},
    "Cu XLPE 750 kcmil 35 kV": {"R1": 0.048, "X1": 0.093, "R0": 0.105, "X0": 0.128, "C1": 0.28},
}

# Standard flat-configuration presets: spacing D (ft) between adjacent phases.
# (Dab = Dbc = D, Dca = 2D for a symmetric horizontal arrangement.)


# ==========================================================================
# SECTION 3 - OVERHEAD TIE LINE (modified Carson)
# ==========================================================================
def carson_line(rac, gmr_ft, OD_in, Dab_ft, Dbc_ft, Dca_ft, height_ft,
                rho, f, length_mi, kv, mva_base=100.0):
    """Return dict of positive/zero sequence R,X (ohm) and B (S), plus per unit."""
    Rd = math.pi**2 * f * 1e-7 * MI                 # earth-resistance term (ohm/mile)
    k = f * MU0 * MI                                # reactance scale (ohm/mile)
    De_ft = (658.37 * math.sqrt(rho / f)) / FT      # Carson depth in ft
    Ke = math.log(De_ft)

    def zself(g):  return complex(rac + Rd, k * (math.log(1.0 / g) + Ke))
    def zmut(D):   return complex(Rd,       k * (math.log(1.0 / D) + Ke))

    zaa = zself(gmr_ft)
    Zm = (zmut(Dab_ft) + zmut(Dbc_ft) + zmut(Dca_ft)) / 3.0
    Z1 = zaa - Zm
    Z0 = zaa + 2.0 * Zm

    # ---- shunt via Maxwell potential coefficients (with earth images) ----
    r_ft = (OD_in / 2.0) / 12.0
    def img(D): return math.sqrt(D**2 + (2.0 * height_ft)**2)
    p_aa = math.log(2.0 * height_ft / r_ft)
    p_ab = math.log(img(Dab_ft) / Dab_ft)
    p_bc = math.log(img(Dbc_ft) / Dbc_ft)
    p_ca = math.log(img(Dca_ft) / Dca_ft)
    ps = p_aa
    pm = (p_ab + p_bc + p_ca) / 3.0
    p1, p0 = ps - pm, ps + 2.0 * pm
    C1 = 2 * math.pi * EPS0 / p1 * MI               # F/mile
    C0 = 2 * math.pi * EPS0 / p0 * MI
    B1 = 2 * math.pi * f * C1                        # S/mile
    B0 = 2 * math.pi * f * C0

    L = length_mi
    R1, X1 = Z1.real * L, Z1.imag * L
    R0, X0 = Z0.real * L, Z0.imag * L
    C1t, C0t = C1 * L, C0 * L
    B1t, B0t = B1 * L, B0 * L

    Zbase = kv**2 / mva_base
    return {
        "R1": R1, "X1": X1, "R0": R0, "X0": X0,
        "C1_uF": C1t * 1e6, "C0_uF": C0t * 1e6, "B1_S": B1t, "B0_S": B0t,
        "R1_pu": R1 / Zbase, "X1_pu": X1 / Zbase,
        "R0_pu": R0 / Zbase, "X0_pu": X0 / Zbase,
        "B1_pu": B1t * Zbase, "B0_pu": B0t * Zbase,
        "Zbase": Zbase, "De_ft": De_ft, "GMD_ft": (Dab_ft*Dbc_ft*Dca_ft)**(1/3.0),
    }


# ==========================================================================
# SECTION 5 - MV COLLECTOR EQUIVALENT (WECC aggregation)
# ==========================================================================
def collector_equiv(S_T, pct_Z, XR, n_per_feeder, n_feeders, kv_c,
                    cable, cable_len_km, x0x1_tx=1.0, f=60.0, mva_base=100.0):
    """
    S_T      transformer rating (MVA)      pct_Z leakage impedance (%)
    XR       transformer X/R ratio         n_per_feeder units per feeder (N)
    n_feeders feeders (M)                  kv_c collector voltage (kV)
    cable    dict R1,X1,R0,X0 (ohm/km), C1 (uF/km)   cable_len_km per-feeder routed length
    x0x1_tx  transformer X0/X1 (1.0 for delta / grounded-wye step-up)
    """
    N, M = n_per_feeder, n_feeders
    # ---- single transformer referred to collector kV ----
    Zb_T = kv_c**2 / S_T
    Zmag = (pct_Z / 100.0) * Zb_T
    X = Zmag / math.sqrt(1.0 + (1.0 / XR)**2)
    R = X / XR
    # ---- aggregate all M*N transformers in parallel ----
    R1_tx = R / (M * N)
    X1_tx = X / (M * N)
    R0_tx = R1_tx                       # Dyn: delta closes zero-seq loop -> Z0~Z1
    X0_tx = X1_tx * x0x1_tx
    # ---- distributed cable: WECC daisy factor, then parallel feeders ----
    fac = (N + 1) * (2 * N + 1) / (6.0 * N**2)
    def cbl_eq(perkm):
        return perkm * cable_len_km * fac / M       # collector-referred, ohm
    R1_c = cbl_eq(cable["R1"]); X1_c = cbl_eq(cable["X1"])
    R0_c = cbl_eq(cable["R0"]); X0_c = cbl_eq(cable["X0"])
    # ---- totals ----
    R1, X1 = R1_tx + R1_c, X1_tx + X1_c
    R0, X0 = R0_tx + R0_c, X0_tx + X0_c
    # ---- shunt (all feeders): C0 = C1 for shielded cable ----
    C1_tot = cable["C1"] * 1e-6 * cable_len_km * M   # F
    B1 = 2 * math.pi * f * C1_tot
    B0 = B1
    Zbase = kv_c**2 / mva_base
    return {
        "R1": R1, "X1": X1, "R0": R0, "X0": X0,
        "C1_uF": C1_tot * 1e6, "C0_uF": C1_tot * 1e6, "B1_S": B1, "B0_S": B0,
        "R1_pu": R1 / Zbase, "X1_pu": X1 / Zbase,
        "R0_pu": R0 / Zbase, "X0_pu": X0 / Zbase,
        "B1_pu": B1 * Zbase, "B0_pu": B0 * Zbase,
        "Zbase": Zbase, "daisy_factor": fac,
        "Z_tx_single_ohm": Zmag, "n_transformers": M * N,
    }
