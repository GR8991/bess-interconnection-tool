import pandas as pd
import numpy as np

# =======================================================================
# BASIC CONDUCTOR DATABASE (COMMON ACSR & AAAC)
# =======================================================================

def load_basic_db():
    """
    Returns a DataFrame containing commonly used transmission & distribution
    conductors with their AC resistance, GMR, and diameter.
    Values are representative typical values used in planning studies.
    Units:
      - R_AC in ohm/km
      - GMR in meters
      - Diameter in meters
    """

    data = {
        "Name": [
            "4/0 ACSR", "336 ACSR Linnet", "477 ACSR Hawk",
            "795 ACSR Drake", "954 ACSR Rail", "1/0 AAAC", "250 AAAC"
        ],
        "R_AC_ohm_per_km": [
            0.272, 0.161, 0.117, 0.082, 0.068, 0.36, 0.23
        ],
        "GMR_m": [
            0.0045, 0.0067, 0.0081, 0.0105, 0.0118, 0.0035, 0.0048
        ],
        "Diameter_m": [
            0.0112, 0.0143, 0.0181, 0.0281, 0.0315, 0.0102, 0.0145
        ],
        "Weight_kg_per_km": [
            430, 640, 820, 1380, 1670, 310, 450
        ]
    }

    df = pd.DataFrame(data)
    return df


# =======================================================================
# FULL IEEE / CIGRE CONDUCTOR LIBRARY
# =======================================================================

def load_full_db():
    """
    Loads the full conductor database from CSV.
    The file must be located at data/conductor_library_full.csv
    If file not found, generates sample dataset.
    """

    try:
        df = pd.read_csv("data/conductor_library_full.csv")
        return df

    except FileNotFoundError:
        # Generate a sample extended database as fallback
        data = {
            "Name": [
                "1350 AAC", "266.8 ACSR Partridge", "397.5 ACSR Ibis",
                "605 ACSR Ostrich", "636 ACSR Grosbeak",
                "795 ACSR Drake", "1033 ACSR Curlew",
                "1113 ACSR Baldfaced", "1272 ACSR Bittern"
            ],
            "R_AC_ohm_per_km": [
                0.292, 0.211, 0.157, 0.107, 0.097, 0.082, 0.067, 0.062, 0.058
            ],
            "GMR_m": [
                0.0032, 0.0073, 0.0089, 0.0107, 0.0119, 0.0105, 0.0132, 0.0141, 0.0151
            ],
            "Diameter_m": [
                0.0105, 0.0175, 0.0201, 0.0245, 0.0279, 0.0281, 0.0317, 0.0331, 0.0348
            ],
            "Weight_kg_per_km": [
                390, 590, 740, 1120, 1300, 1380, 1540, 1620, 1710
            ]
        }

        df = pd.DataFrame(data)
        return df


# =======================================================================
# UTILITY FUNCTIONS
# =======================================================================

def get_conductor_parameters(name):
    """
    Retrieve conductor parameters from either basic or full DB.
    Returns dict with r_ac, gmr, diameter.
    """

    df_basic = load_basic_db()
    df_full = load_full_db()

    if name in df_basic["Name"].values:
        row = df_basic[df_basic["Name"] == name].iloc[0]
    elif name in df_full["Name"].values:
        row = df_full[df_full["Name"] == name].iloc[0]
    else:
        raise ValueError(f"Conductor '{name}' not found in database.")

    return {
        "r_ac": row["R_AC_ohm_per_km"],
        "gmr": row["GMR_m"],
        "diameter": row["Diameter_m"]
    }
