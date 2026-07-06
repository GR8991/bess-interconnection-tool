"""
Interconnection Impedance Tool  -  Section 3 (Tie Line) & Section 5 (Collector)
IEEE/ANSI methods: modified Carson's equations + WECC collector aggregation.

Run:   pip install streamlit
       streamlit run app.py
"""
import streamlit as st
import ieee_calc as ic

st.set_page_config(page_title="IBR Interconnection Impedance Tool", layout="wide")

M_PER_FT = 0.3048

def to_ft(value, unit):
    return value / M_PER_FT if unit == "m" else value

def len_to_mi(value, unit):
    return {"mi": value, "km": value / 1.60934, "ft": value / 5280.0}[unit]

def len_to_km(value, unit):
    return {"km": value, "mi": value * 1.60934, "ft": value * 0.0003048}[unit]

st.title("IBR Interconnection Impedance Tool")
st.caption("Section 3 - Interconnection Tie Line  |  Section 5 - Collector System Equivalent   "
           "(modified Carson / WECC aggregation, IEEE-ANSI practice)")

tab3, tab5, about = st.tabs(["Section 3 - Tie Line",
                             "Section 5 - Collector System",
                             "About / Method"])

# =====================================================================
# SECTION 3
# =====================================================================
with tab3:
    st.subheader("Overhead transmission tie line")
    c1, c2, c3 = st.columns(3)
    with c1:
        cond_name = st.selectbox("Conductor (ACSR)", list(ic.CONDUCTORS.keys()),
                                 index=6, key="cond")
        cond = ic.CONDUCTORS[cond_name]
        kv = st.number_input("Line voltage (kV)", 1.0, 800.0, 100.0, 1.0, key="kv3")
        length_val = st.number_input("Line length", 0.0, 10000.0, 1.0, 0.1, key="len3")
        length_unit = st.radio("Length unit", ["mi", "km", "ft"], horizontal=True, key="lu3")
    with c2:
        geom_unit = st.radio("Spacing / height unit", ["m", "ft"], horizontal=True, key="gu3")
        D = st.number_input(f"Adjacent phase spacing D ({geom_unit})",
                            0.1, 100.0, 7.0, 0.1, key="D3")
        height = st.number_input(f"Conductor height ({geom_unit})",
                                 1.0, 200.0, 12.0, 0.5, key="h3")
        config = st.selectbox("Configuration",
                              ["Flat horizontal (Dca = 2D)", "Custom spacings"], key="cfg3")
    with c3:
        rho = st.number_input("Earth resistivity (ohm-m)", 1.0, 10000.0, 100.0, 10.0, key="rho3")
        f = st.number_input("Frequency (Hz)", 50.0, 60.0, 60.0, 10.0, key="f3")
        mva = st.number_input("MVA base", 1.0, 100000.0, 100.0, 10.0, key="mva3")

    with st.expander("Override conductor data / custom spacings"):
        o1, o2, o3 = st.columns(3)
        rac = o1.number_input("AC resistance (ohm/mile)", 0.0, 5.0, float(cond["rac"]),
                              0.0001, format="%.4f", key="rac3")
        gmr = o2.number_input("GMR (ft)", 0.0001, 1.0, float(cond["gmr_ft"]),
                              0.0001, format="%.4f", key="gmr3")
        od = o3.number_input("Outside diameter (in)", 0.01, 5.0, float(cond["OD_in"]),
                             0.001, format="%.3f", key="od3")
        if config == "Custom spacings":
            s1, s2, s3 = st.columns(3)
            Dab = to_ft(s1.number_input(f"Dab ({geom_unit})", 0.1, 100.0, D, 0.1, key="Dab"), geom_unit)
            Dbc = to_ft(s2.number_input(f"Dbc ({geom_unit})", 0.1, 100.0, D, 0.1, key="Dbc"), geom_unit)
            Dca = to_ft(s3.number_input(f"Dca ({geom_unit})", 0.1, 100.0, 2 * D, 0.1, key="Dca"), geom_unit)
        else:
            Dab = Dbc = to_ft(D, geom_unit); Dca = to_ft(2 * D, geom_unit)

    if config == "Flat horizontal (Dca = 2D)":
        Dab = Dbc = to_ft(D, geom_unit); Dca = to_ft(2 * D, geom_unit)

    res = ic.carson_line(rac=rac, gmr_ft=gmr, OD_in=od,
                         Dab_ft=Dab, Dbc_ft=Dbc, Dca_ft=Dca,
                         height_ft=to_ft(height, geom_unit),
                         rho=rho, f=f, length_mi=len_to_mi(length_val, length_unit),
                         kv=kv, mva_base=mva)

    st.markdown("#### Results - Section 3 (c)-(h)")
    rows = [
        ("(c) R1  positive-seq resistance", f"{res['R1']:.4f} ohm", f"{res['R1_pu']:.6f} pu"),
        ("(d) X1  positive-seq reactance",  f"{res['X1']:.4f} ohm", f"{res['X1_pu']:.6f} pu"),
        ("(e) B1 / C1  positive-seq shunt", f"C1 = {res['C1_uF']:.4f} uF", f"B1 = {res['B1_pu']:.3e} pu"),
        ("(f) R0  zero-seq resistance",     f"{res['R0']:.4f} ohm", f"{res['R0_pu']:.6f} pu"),
        ("(g) X0  zero-seq reactance",      f"{res['X0']:.4f} ohm", f"{res['X0_pu']:.6f} pu"),
        ("(h) B0 / C0  zero-seq shunt",     f"C0 = {res['C0_uF']:.4f} uF", f"B0 = {res['B0_pu']:.3e} pu"),
    ]
    st.table({"Form item": [r[0] for r in rows],
              "Ohms / F":  [r[1] for r in rows],
              "Per unit (on base)": [r[2] for r in rows]})
    st.caption(f"Z_base = {res['Zbase']:.2f} ohm  |  Carson depth De = {res['De_ft']:.0f} ft  "
               f"|  GMD = {res['GMD_ft']:.2f} ft. Values are for the entered length.")

    txt3 = ("Section 3 - Tie Line\n"
            f"Conductor: {cond_name}\nLine voltage: {kv} kV   Length: {length_val} {length_unit}\n"
            f"R1={res['R1']:.4f} ohm ({res['R1_pu']:.6f} pu)\n"
            f"X1={res['X1']:.4f} ohm ({res['X1_pu']:.6f} pu)\n"
            f"C1={res['C1_uF']:.4f} uF  B1={res['B1_pu']:.3e} pu\n"
            f"R0={res['R0']:.4f} ohm ({res['R0_pu']:.6f} pu)\n"
            f"X0={res['X0']:.4f} ohm ({res['X0_pu']:.6f} pu)\n"
            f"C0={res['C0_uF']:.4f} uF  B0={res['B0_pu']:.3e} pu\n")
    st.download_button("Download Section 3 values (.txt)", txt3, "section3_tieline.txt", key="dl3")

# =====================================================================
# SECTION 5
# =====================================================================
with tab5:
    st.subheader("MV collector system equivalent")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**Inverter step-up transformer**")
        S_T = st.number_input("Rating (MVA each)", 0.1, 100.0, 5.0, 0.1, key="st5")
        pctZ = st.number_input("Leakage impedance (%)", 1.0, 20.0, 7.0, 0.1, key="pz5")
        XR = st.number_input("X/R ratio", 1.0, 50.0, 10.0, 0.5, key="xr5")
        conn = st.selectbox("Winding connection",
                            ["Delta LV / Wye-g HV  (X0=X1)",
                             "Wye-g / Wye-g  (enter X0/X1)",
                             "Custom X0/X1"], key="conn5")
    with c2:
        st.markdown("**Layout**")
        N = st.number_input("Transformers per feeder (N)", 1, 50, 5, 1, key="n5")
        Mf = st.number_input("Number of feeders (M)", 1, 50, 6, 1, key="m5")
        kv_c = st.number_input("Collector voltage (kV)", 1.0, 69.0, 34.5, 0.5, key="kvc5")
        mva5 = st.number_input("MVA base", 1.0, 100000.0, 100.0, 10.0, key="mva5")
    with c3:
        st.markdown("**MV cable (per feeder)**")
        cab_name = st.selectbox("Cable", list(ic.CABLES.keys()), index=5, key="cab5")
        cab = dict(ic.CABLES[cab_name])
        clen = st.number_input("Routed feeder length", 0.0, 100000.0, 1500.0, 10.0, key="cl5")
        clen_unit = st.radio("Cable length unit", ["ft", "m", "km"], horizontal=True, key="clu5")

    x0x1 = 1.0
    if conn != "Delta LV / Wye-g HV  (X0=X1)":
        x0x1 = st.slider("Transformer X0 / X1", 0.5, 3.0, 1.0, 0.05, key="x0x15")

    with st.expander("Override cable per-length data (ohm/km, uF/km)"):
        o = st.columns(5)
        cab["R1"] = o[0].number_input("R1", 0.0, 5.0, float(cab["R1"]), 0.001, format="%.3f", key="cR1")
        cab["X1"] = o[1].number_input("X1", 0.0, 5.0, float(cab["X1"]), 0.001, format="%.3f", key="cX1")
        cab["R0"] = o[2].number_input("R0", 0.0, 5.0, float(cab["R0"]), 0.001, format="%.3f", key="cR0")
        cab["X0"] = o[3].number_input("X0", 0.0, 5.0, float(cab["X0"]), 0.001, format="%.3f", key="cX0")
        cab["C1"] = o[4].number_input("C1", 0.0, 2.0, float(cab["C1"]), 0.001, format="%.3f", key="cC1")

    clen_km = len_to_km(clen, clen_unit)
    res5 = ic.collector_equiv(S_T=S_T, pct_Z=pctZ, XR=XR, n_per_feeder=int(N),
                              n_feeders=int(Mf), kv_c=kv_c, cable=cab,
                              cable_len_km=clen_km, x0x1_tx=x0x1, mva_base=mva5)

    st.markdown("#### Results - Section 5 (c) i-vi")
    rows5 = [
        ("(i) R1  positive-seq resistance", f"{res5['R1']:.4f} ohm", f"{res5['R1_pu']:.6f} pu"),
        ("(ii) X1  positive-seq reactance", f"{res5['X1']:.4f} ohm", f"{res5['X1_pu']:.6f} pu"),
        ("(iii) B1 / C1 positive-seq shunt", f"C1 = {res5['C1_uF']:.4f} uF", f"B1 = {res5['B1_pu']:.3e} pu"),
        ("(iv) R0  zero-seq resistance",    f"{res5['R0']:.4f} ohm", f"{res5['R0_pu']:.6f} pu"),
        ("(v) X0  zero-seq reactance",      f"{res5['X0']:.4f} ohm", f"{res5['X0_pu']:.6f} pu"),
        ("(vi) B0 / C0 zero-seq shunt",     f"C0 = {res5['C0_uF']:.4f} uF", f"B0 = {res5['B0_pu']:.3e} pu"),
    ]
    st.table({"Form item": [r[0] for r in rows5],
              "Ohms / F":  [r[1] for r in rows5],
              "Per unit (on base)": [r[2] for r in rows5]})
    st.caption(f"Z_base = {res5['Zbase']:.2f} ohm  |  {res5['n_transformers']} transformers aggregated  "
               f"|  WECC daisy factor = {res5['daisy_factor']:.3f}  "
               f"|  single-transformer |Z| = {res5['Z_tx_single_ohm']:.2f} ohm.")

    txt5 = ("Section 5 - Collector System Equivalent\n"
            f"Transformer: {S_T} MVA {pctZ}% X/R={XR}  x{res5['n_transformers']} units\n"
            f"Cable: {cab_name}   feeder length: {clen} {clen_unit}\n"
            f"Collector voltage: {kv_c} kV\n"
            f"R1={res5['R1']:.4f} ohm ({res5['R1_pu']:.6f} pu)\n"
            f"X1={res5['X1']:.4f} ohm ({res5['X1_pu']:.6f} pu)\n"
            f"C1={res5['C1_uF']:.4f} uF  B1={res5['B1_pu']:.3e} pu\n"
            f"R0={res5['R0']:.4f} ohm ({res5['R0_pu']:.6f} pu)\n"
            f"X0={res5['X0']:.4f} ohm ({res5['X0_pu']:.6f} pu)\n"
            f"C0={res5['C0_uF']:.4f} uF  B0={res5['B0_pu']:.3e} pu\n")
    st.download_button("Download Section 5 values (.txt)", txt5, "section5_collector.txt", key="dl5")

# =====================================================================
# ABOUT
# =====================================================================
with about:
    st.markdown("""
### Method
**Section 3 - Tie line (modified Carson's equations).**
Self `z_ii = (r_ac + R_d) + jk[ln(1/GMR)+K_e]` and mutual `z_ij = R_d + jk[ln(1/D_ij)+K_e]`,
with `R_d = pi^2 f 1e-7 * 1609`, `k = f mu0 * 1609`, `K_e = ln(D_e)`, `D_e = 658.4 sqrt(rho/f)` m.
Sequence: `Z1 = Zs - Zm`, `Z0 = Zs + 2 Zm`. Shunt from Maxwell potential-coefficient matrix
(earth images); `B = wC`.

**Section 5 - Collector (WECC aggregation).**
Transformer `|Z| = (%Z/100) * kV^2/MVA`, split by X/R, then all `M*N` units in parallel.
Distributed feeder cable uses the daisy-chain factor `(N+1)(2N+1)/(6N^2)`; feeders combine in
parallel. For a delta / grounded-wye step-up `Z0 = Z1`. For shielded cable `C0 = C1`.

Per unit: `series pu = ohm / Zbase`, `shunt pu = B * Zbase`, `Zbase = kV^2 / MVA_base`.

### Notes
- Library resistance values are typical (75 C); override with the datasheet.
- Cable per-length values are representative; override with IEC 60287 / manufacturer data.
- Use the actual routed feeder length - shunt susceptance scales directly with length.
""")
