import numpy as np
import pandas as pd


# =====================================================================
# CONSTANTS
# =====================================================================
MU0 = 4 * np.pi * 1e-7               # H/m
EPS0 = 8.854e-12                     # F/m
FREQ = 60                            # Hz
OMEGA = 2 * np.pi * FREQ


# =====================================================================
# IMPEDANCE BASE CALCULATIONS
# =====================================================================
def base_impedance(kV_base, MVA_base):
    """
    Compute impedance base Z_base = V^2 / S.
    kV_base: line-to-line base (kV)
    MVA_base: 3-phase base power
    Returns Z_base in ohms.
    """
    return (kV_base ** 2) / MVA_base


def pu_conversion(Z_ohm, Z_base):
    """
    Convert impedance in ohms to per-unit.
    """
    if Z_base == 0:
        return 0
    return Z_ohm / Z_base


def ohm_conversion(Z_pu, Z_base):
    """
    Convert per-unit impedance to ohms.
    """
    return Z_pu * Z_base


# =====================================================================
# LENGTH UNIT CONVERSIONS
# =====================================================================
def km_to_miles(km):
    return km * 0.621371

def miles_to_km(miles):
    return miles * 1.60934

def convert_length_to_meters(length, unit):
    if unit.lower() == "km":
        return length * 1000
    elif unit.lower() == "mile":
        return length * 1609.34
    else:
        return length


# =====================================================================
# JSON FLATTENER FOR EXPORT
# =====================================================================
def flatten_dict(d, parent_key='', sep='_'):
    """
    Flattens nested dictionaries into a flat structure for DataFrame export.
    """
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k

        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))

    return dict(items)


# =====================================================================
# SAFE NUMBER FORMATTER
# =====================================================================
def fmt(x, decimals=3):
    """
    Safely format a number for display.
    """
    try:
        return round(float(x), decimals)
    except:
        return x


# =====================================================================
# USED BY LINE MODELS
# =====================================================================
def impedance_to_mag_angle(Z):
    """
    Convert complex impedance to magnitude and angle.
    Returns (mag, angle_deg)
    """
    return np.abs(Z), np.angle(Z, deg=True)
