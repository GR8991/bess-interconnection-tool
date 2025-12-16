import streamlit as st
import pandas as pd
import numpy as np

from modules.line_calc import compute_overhead_impedance, compute_underground_impedance, compute_multicircuit
from modules.conductor_db import load_basic_db, load_full_db
from modules.interconnection_form import interconnection_summary
from modules.export_utils import export_csv, export_excel, export_pdf
from modules.utils import pu_conversion, base_impedance

# --------------------------------------------------------
# ----------- Streamlit Page Config + Minimal CSS --------
# --------------------------------------------------------
st.set_page_config(
    page_title="BESS Interconnection Engineering Tool",
    layout="wide"
)

# Minimal hybrid CSS styling
st.markdown("""
<style>
    .main {background-color: #F8F9FB;}
    .stTabs [role="tablist"] button {
        font-size: 16px;
        padding: 8px 16px;
    }
    .card {
        background-color: white;
        padding: 18px;
        border-radius: 8px;
        box-shadow: 0px 1px 4px rgba(0,0,0,0.1);
        margin-bottom: 15px;
    }
</style>
""", unsafe_allow_html=True)

# --------------------------------------------------------
# -------------------- Sidebar Menu ----------------------
# --------------------------------------------------------
tabs = [
    "Overhead Line Calculator",
    "Underground Cable Calculator",
    "Multi-Circuit Line Model",
    "Per-Unit Conversion",
    "Conductor Database",
    "Interconnection Form Auto-Fill",
    "Export Results"
]

page = st.sidebar.radio("Navigation", tabs)



# =====================================================================
# ====================== TAB 1: OVERHEAD LINE CALCULATOR ==============
# =====================================================================
if page == "Overhead Line Calculator":

    st.title("Overhead Line Impedance Calculator (Carson's Equations)")

    with st.container():
        st.subheader("Input Parameters")

        col1, col2, col3 = st.columns(3)

        conductor_type = col1.selectbox("Select Conductor Type", load_basic_db()["Name"].tolist())
        length_value = col2.number_input("Line Length", value=1.0, min_value=0.0)
        length_unit = col3.selectbox("Length Unit", ["km", "mile"])

        Dab = col1.number_input("Phase Spacing D_ab (m)", value=6.0)
        Dbc = col2.number_input("Phase Spacing D_bc (m)", value=6.0)
        Dca = col3.number_input("Phase Spacing D_ca (m)", value=12.0)

        ha = col1.number_input("Height of Phase A (m)", value=12.0)
        hb = col2.number_input("Height of Phase B (m)", value=12.0)
        hc = col3.number_input("Height of Phase C (m)", value=12.0)

        rho = st.number_input("Earth Resistivity (Ω·m)", value=100.0)

        if st.button("Compute Impedance"):
            df = load_basic_db()
            row = df[df["Name"] == conductor_type].iloc[0]

            Z1, Z0, B1, B0 = compute_overhead_impedance(
                r_ac=row["R_AC_ohm_per_km"],
                gmr=row["GMR_m"],
                diameter=row["Diameter_m"],
                Dab=Dab, Dbc=Dbc, Dca=Dca,
                ha=ha, hb=hb, hc=hc,
                rho=rho,
                length=length_value,
                unit=length_unit
            )

            st.subheader("Results (Ohms & Per-Unit)")

            st.write("### Positive Sequence (Z1)")
            st.json(Z1)

            st.write("### Zero Sequence (Z0)")
            st.json(Z0)

            st.write("### Shunt Susceptance (B1, B0)")
            st.json({"B1": B1, "B0": B0})



# =====================================================================
# ===================== TAB 2: UNDERGROUND CABLE MODEL ================
# =====================================================================
elif page == "Underground Cable Calculator":

    st.title("Underground Cable Impedance Calculator")

    with st.container():
        st.subheader("Cable Inputs")

        col1, col2, col3 = st.columns(3)

        cond_r = col1.number_input("Conductor Resistance R (Ω/km)", value=0.08)
        cond_gmr = col2.number_input("GMR (m)", value=0.01)
        cond_diam = col3.number_input("Conductor Diameter (m)", value=0.02)

        spacing = col1.number_input("Cable Spacing (m)", value=0.1)
        sheath_bond = col2.selectbox("Sheath Bonding", ["Both Ends Bonded", "Single-Point", "Cross-Bonded"])

        length = col3.number_input("Cable Length (km)", value=1.0)

        if st.button("Compute Cable Impedance"):
            Z1, Z0, B1, B0 = compute_underground_impedance(
                cond_r, cond_gmr, cond_diam, spacing, sheath_bond, length
            )

            st.subheader("Cable Impedance Output")
            st.json({
                "Positive Sequence Z1": Z1,
                "Zero Sequence Z0": Z0,
                "B1": B1,
                "B0": B0
            })



# =====================================================================
# ===================== TAB 3: MULTI-CIRCUIT LINE MODEL ===============
# =====================================================================
elif page == "Multi-Circuit Line Model":

    st.title("Multi-Circuit Overhead Line Modeling")

    st.info("Supports double-circuit, parallel lines, and ground wire modeling.")

    with st.container():
        col1, col2, col3 = st.columns(3)

        circuits = col1.number_input("Number of Circuits", value=2, min_value=1)
        ground_wire = col2.checkbox("Include Overhead Ground Wire?")
        frequency = col3.number_input("Frequency (Hz)", value=60.0)

        if st.button("Compute Multi-Circuit Impedance"):
            results = compute_multicircuit(circuits, ground_wire, frequency)
            st.json(results)



# =====================================================================
# ========================= TAB 4: PER-UNIT CONVERSION ================
# =====================================================================
elif page == "Per-Unit Conversion":

    st.title("Per-Unit Conversion Tool")

    col1, col2, col3 = st.columns(3)

    Z_mag = col1.number_input("Impedance Magnitude (Ω)", value=1.0)
    V_kV = col2.number_input("Base Voltage (kV)", value=34.5)
    MVA = col3.number_input("Base MVA", value=100.0)

    if st.button("Convert to Per-Unit"):
        Z_base = base_impedance(V_kV, MVA)
        Z_pu = pu_conversion(Z_mag, Z_base)
        st.success(f"Per-Unit Impedance = {Z_pu:.6f} pu")



# =====================================================================
# ========================= TAB 5: CONDUCTOR DATABASE =================
# =====================================================================
elif page == "Conductor Database":

    st.title("Conductor Database (Basic + IEEE Full Library)")

    option = st.selectbox("Choose Library", ["Basic", "Full IEEE"])

    if option == "Basic":
        df = load_basic_db()
    else:
        df = load_full_db()

    st.dataframe(df)



# =====================================================================
# ================ TAB 6: INTERCONNECTION FORM AUTO-FILL ==============
# =====================================================================
elif page == "Interconnection Form Auto-Fill":

    st.title("Automatic Interconnection Form Generator")

    st.info("Fill in the BESS project data below to auto-generate a complete interconnection summary.")

    result = interconnection_summary()

    if result:
        st.subheader("Interconnection Form Summary")
        st.json(result)



# =====================================================================
# =========================== TAB 7: EXPORT RESULTS ===================
# =====================================================================
elif page == "Export Results":

    st.title("Export Results")

    st.warning("Generate impedance or interconnection results in other tabs first.")

    df = st.session_state.get("export_df", None)
    if df is None:
        st.error("No data available for export.")
    else:
        if st.button("Export CSV"):
            export_csv(df)

        if st.button("Export Excel"):
            export_excel(df)

        if st.button("Export PDF"):
            export_pdf(df)


