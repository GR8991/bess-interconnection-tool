import streamlit as st
import numpy as np
from modules.utils import base_impedance, pu_conversion

# =====================================================================
#  Interconnection Form Auto-Fill Engine
# =====================================================================

def interconnection_summary():
    """
    Main function that renders the input UI in Streamlit
    and returns a dictionary containing all formatted results.
    """

    st.header("BESS Interconnection Form Auto-Fill")

    # ------------------------------------------------------------
    # ---------------- Input Panel -------------------------------
    # ------------------------------------------------------------
    st.subheader("ESR (BESS) Parameters")

    col1, col2, col3 = st.columns(3)
    poimw = col1.number_input("POI MW Rating", value=100.0)
    poimwh = col2.number_input("POI MWh Rating", value=200.0)
    pcs_count = col3.number_input("Number of PCS/Inverters", value=60)

    pcs_kw = col1.number_input("PCS Rating (kW per PCS)", value=2500.0)
    pcs_kva = col2.number_input("PCS Rating (kVA per PCS)", value=2520.0)
    batt_mwh_gross = col3.number_input("Gross Installed Energy (kWh)", value=294000.0)

    st.subheader("GSU Transformer Parameters")

    colA, colB, colC = st.columns(3)
    gsu_mva_onan = colA.number_input("GSU Rating (ONAN MVA)", value=125.0)
    gsu_mva_onaf = colB.number_input("GSU Rating (ONAF MVA)", value=166.0)
    gsu_impedance = colC.number_input("GSU Positive Sequence Impedance (%)", value=10.0)

    # Transformer voltages
    hv_kv = colA.number_input("GSU HV Voltage (kV)", value=115.0)
    mv_kv = colB.number_input("GSU MV Voltage (kV)", value=34.5)
    tv_kv = colC.number_input("GSU Tertiary Voltage (kV)", value=13.8)

    # Tertiary MVA
    tertiary_mva = colA.number_input("Tertiary Rating (MVA)", value=25.0)

    # Tap positions
    st.subheader("GSU Tap Positions (kV)")
    taps = compute_tap_positions(hv_kv)

    # ------------------------------------------------------------
    # Collector System Equivalent
    # ------------------------------------------------------------
    st.subheader("Collector System Equivalent Model (34.5 kV Bus)")

    colD, colE, colF = st.columns(3)
    coll_mva = colD.number_input("Collector System Rating (MVA)", value=125.0)
    coll_r1 = colE.number_input("R1 (pu)", value=0.002)
    coll_x1 = colF.number_input("X1 (pu)", value=0.06)
    coll_b1 = colD.number_input("B1 (pu)", value=0.30)
    coll_r0 = colE.number_input("R0 (pu)", value=0.003)
    coll_x0 = colF.number_input("X0 (pu)", value=0.10)
    coll_b0 = colD.number_input("B0 (pu)", value=0.30)

    # ------------------------------------------------------------
    # Reactive Power Capability
    # ------------------------------------------------------------
    st.subheader("Reactive Power Capability")

    q_per_pcs = colA.number_input("Reactive Capability per PCS (MVAR)", value=0.56)

    total_q = q_per_pcs * pcs_count

    # ------------------------------------------------------------
    # Compute ESR Gross & Charging Capability
    # ------------------------------------------------------------
    esr_gross_kw = pcs_kw * pcs_count
    esr_gross_kva = pcs_kva * pcs_count

    summary = {
        "ESR Gross Capability": {
            "Peak Power (kW)": esr_gross_kw,
            "Max Energy (kWh)": batt_mwh_gross
        },
        "ESR Charging Capability": {
            "Peak Charging Power (kW)": esr_gross_kw,
            "Max Energy (kWh)": batt_mwh_gross
        },
        "PCS Fleet": {
            "PCS Count": pcs_count,
            "PCS kW": pcs_kw,
            "PCS kVA": pcs_kva
        },
        "Reactive Capability": {
            "Per PCS (MVAR)": q_per_pcs,
            "Total Plant Reactive Capability (MVAR)": total_q
        },
        "GSU Transformer": {
            "Rating ONAN (MVA)": gsu_mva_onan,
            "Rating ONAF (MVA)": gsu_mva_onaf,
            "HV (kV)": hv_kv,
            "MV (kV)": mv_kv,
            "Tertiary (kV)": tv_kv,
            "Tertiary Rating (MVA)": tertiary_mva,
            "Impedance (%)": gsu_impedance,
            "Tap Positions (kV)": taps,
            "BIL Ratings (kV)": {
                "HV": 550,
                "MV": 230,
                "Tertiary": 110
            }
        },
        "Collector System Equivalent": {
            "Voltage (kV)": 34.5,
            "Rating (MVA)": coll_mva,
            "R1 (pu)": coll_r1,
            "X1 (pu)": coll_x1,
            "B1 (pu)": coll_b1,
            "R0 (pu)": coll_r0,
            "X0 (pu)": coll_x0,
            "B0 (pu)": coll_b0
        },
        "POI Data": {
            "POI MW": poimw,
            "POI MWh": poimwh,
            "POI Voltage (kV)": hv_kv,
            "Required PF": 0.95
        }
    }

    return summary


# =====================================================================
#  Helper: Compute GSU Tap Positions in kV
# =====================================================================
def compute_tap_positions(hv_kv):
    """
    Compute tap positions at -5%, -2.5%, 0%, +2.5%, +5%
    """
    tap_minus_5 = hv_kv * 0.95
    tap_minus_2_5 = hv_kv * 0.975
    tap_nominal = hv_kv
    tap_plus_2_5 = hv_kv * 1.025
    tap_plus_5 = hv_kv * 1.05

    return {
        "-5%": round(tap_minus_5, 3),
        "-2.5%": round(tap_minus_2_5, 3),
        "0% (Nominal)": round(tap_nominal, 3),
        "+2.5%": round(tap_plus_2_5, 3),
        "+5%": round(tap_plus_5, 3)
    }
