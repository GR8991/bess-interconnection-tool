def compute_overhead_impedance(
    r_ac_km, gmr_m, diameter_m,
    Dab, Dbc, Dca,
    ha, hb, hc,
    rho,
    length,
    unit="mile"
):
    """
    Returns sequence impedances for an overhead line using Carson's method.
    r_ac_km: AC resistance in ohm/km
    gmr_m: geometric mean radius (m)
    diameter_m: conductor diameter (m)
    length: line length in mile or km
    """

    # Convert r_ac from ohm/km → ohm/meter
    r_ac = r_ac_km / 1000.0

    # Convert line length to meters
    if unit.lower() == "mile":
        L = length * 1609.34
    else:
        L = length * 1000

    # Earth correction constant
    k = (1 + 1j) * 0.000001 * np.sqrt(OMEGA * MU0 / (2 * rho))

    # Self impedance
    def Zself(h):
        return r_ac + 1j * OMEGA * MU0 / (2 * np.pi) * (np.log(2 * h / gmr_m) + k)

    # Mutual impedance
    def Zmutual(D, h1, h2):
        De = np.sqrt(D * D + (h1 + h2) * (h1 + h2))
        return 1j * OMEGA * MU0 / (2 * np.pi) * (np.log(De / D) + k)

    # Build primitive Z matrix
    Z = np.zeros((3, 3), dtype=complex)

    heights = [ha, hb, hc]

    for i in range(3):
        Z[i, i] = Zself(heights[i])

    Z[0, 1] = Z[1, 0] = Zmutual(Dab, ha, hb)
    Z[1, 2] = Z[2, 1] = Zmutual(Dbc, hb, hc)
    Z[0, 2] = Z[2, 0] = Zmutual(Dca, ha, hc)

    # Symmetrical component transform
    a = np.exp(1j * 2 * np.pi / 3)
    T = (1 / 3) * np.array([
        [1, 1, 1],
        [1, a, a ** 2],
        [1, a ** 2, a]
    ])

    Z012 = T @ Z @ np.linalg.inv(T)

    # Sequence impedances per meter
    Z1_per_m = Z012[1, 1]
    Z0_per_m = Z012[0, 0]

    # Convert to total ohms (multiply by line length)
    Z1 = Z1_per_m * L
    Z0 = Z0_per_m * L

    # Shunt susceptance based on standard overhead C'
    C_per_m = (2 * np.pi * EPS0) / np.log(Dab / (diameter_m / 2))
    B1 = OMEGA * C_per_m * L
    B0 = B1 * 0.75  # zero seq slightly smaller

    Z1_dict = {
        "R1_ohm": np.real(Z1),
        "X1_ohm": np.imag(Z1),
        "Z1_mag": np.abs(Z1),
        "Z1_angle_deg": np.angle(Z1, deg=True),
    }

    Z0_dict = {
        "R0_ohm": np.real(Z0),
        "X0_ohm": np.imag(Z0),
        "Z0_mag": np.abs(Z0),
        "Z0_angle_deg": np.angle(Z0, deg=True),
    }

    return Z1_dict, Z0_dict, B1, B0
