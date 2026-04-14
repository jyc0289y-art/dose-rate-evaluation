import math

def log_interp(x, x1, y1, x2, y2):
    t = (math.log(x) - math.log(x1)) / (math.log(x2) - math.log(x1))
    return math.exp(math.log(y1) + t * (math.log(y2) - math.log(y1)))

# NIST Air mu_en/rho (keV, cm^2/g)
air_data = [
    (5, 39.31), (6, 22.70), (8, 9.446), (10, 4.742),
    (15, 1.334), (20, 0.5389), (30, 0.1537), (40, 0.06833),
    (50, 0.04098), (60, 0.03041), (80, 0.02407), (100, 0.02325)
]

def get_air_muen(E):
    # [Gemini 교차검증 수정] 원래 코드는 E < air_data[0][0] (5 keV 미만) 입력 시
    # air_data[-1][1] (100 keV 값, 0.02325)을 반환하는 fallback 버그가 있었음.
    # 저에너지일수록 μen/ρ이 급격히 증가하므로 100 keV 값 반환은 심각한 과소평가.
    # 수정: 처음 두 점(5, 6 keV)으로 log-log 외삽.
    if E < air_data[0][0]:
        e1, m1 = air_data[0]
        e2, m2 = air_data[1]
        return log_interp(E, e1, m1, e2, m2)
    for i in range(len(air_data)-1):
        e1, m1 = air_data[i]
        e2, m2 = air_data[i+1]
        if e1 <= E <= e2:
            return log_interp(E, e1, m1, e2, m2)
    return air_data[-1][1]

def show_interp_detail(E, label):
    for i in range(len(air_data)-1):
        e1, m1 = air_data[i]
        e2, m2 = air_data[i+1]
        if e1 <= E <= e2:
            t = (math.log(E) - math.log(e1)) / (math.log(e2) - math.log(e1))
            ln_y = math.log(m1) + t * (math.log(m2) - math.log(m1))
            result = math.exp(ln_y)
            print(f"  {label} ({E} keV):")
            print(f"    NIST grid: {e1} keV ({m1}) ~ {e2} keV ({m2})")
            print(f"    t = [ln({E})-ln({e1})] / [ln({e2})-ln({e1})]")
            print(f"      = [{math.log(E):.4f}-{math.log(e1):.4f}] / [{math.log(e2):.4f}-{math.log(e1):.4f}]")
            print(f"      = {math.log(E)-math.log(e1):.4f} / {math.log(e2)-math.log(e1):.4f} = {t:.4f}")
            print(f"    ln(mu) = ln({m1}) + {t:.4f} * [ln({m2})-ln({m1})]")
            print(f"           = {math.log(m1):.4f} + {t:.4f} * ({math.log(m2):.4f}-{math.log(m1):.4f})")
            print(f"           = {math.log(m1):.4f} + {t:.4f} * ({math.log(m2)-math.log(m1):.4f})")
            print(f"           = {ln_y:.4f}")
            print(f"    mu_en/rho = e^{ln_y:.4f} = {result:.4f} cm^2/g")
            return result
    return None

# Am-241 emissions (ICRP 107)
am241_lines = [
    (11.89, 0.0083, "Np Ll"),
    (13.27, 0.0148, "Np La2"),
    (13.95, 0.1303, "Np La1"),
    (16.84, 0.0481, "Np Lb2"),
    (17.06, 0.006, "Np Lb3"),
    (17.75, 0.1896, "Np Lb1"),
    (20.78, 0.0488, "Np Lg1"),
    (21.10, 0.016, "Np Lg2,3"),
    (26.34, 0.024, "g 26.3"),
    (59.54, 0.3592, "g 59.5"),
]

# Cd-109 emissions (ICRP 107)
cd109_lines = [
    (21.99, 0.290, "Ag Ka2"),
    (22.16, 0.553, "Ag Ka1"),
    (24.91, 0.106, "Ag Kb1,3"),
    (25.46, 0.043, "Ag Kb2"),
    (88.03, 0.0363, "g 88.0"),
]

keV_to_J = 1.602e-16  # J/keV
cm2g_to_m2kg = 0.1

for name, lines in [("Am-241", am241_lines), ("Cd-109", cd109_lines)]:
    print("=" * 80)
    print(f"  {name} Gamma Constant Calculation from NIST Data")
    print("=" * 80)

    # Step 1: Show all interpolations
    print(f"\n--- Step 1: Air mu_en/rho interpolation ---")
    muens = {}
    for E, Y, label in lines:
        mu = show_interp_detail(E, label)
        muens[E] = mu
        print()

    # Step 2: Calculate each line's contribution
    print(f"\n--- Step 2: Per-line contribution Y_i * E_i[J] * (mu_en/rho)[m^2/kg] ---")
    total_sum = 0
    contributions = []
    for E, Y, label in lines:
        E_J = E * keV_to_J
        mu_m2kg = muens[E] * cm2g_to_m2kg
        contrib = Y * E_J * mu_m2kg
        total_sum += contrib
        contributions.append((E, Y, label, muens[E], E_J, mu_m2kg, contrib))
        print(f"  {label} ({E} keV):")
        print(f"    Y={Y}, E={E_J:.4e} J, mu_en/rho={mu_m2kg:.4f} m^2/kg")
        print(f"    contrib = {Y} * {E_J:.4e} * {mu_m2kg:.4f} = {contrib:.4e}")

    print(f"\n  Sum = {total_sum:.4e}")

    # Step 3: Divide by 4pi
    gamma_SI = total_sum / (4 * math.pi)
    print(f"\n--- Step 3: Gamma_SI = Sum / 4pi ---")
    print(f"  = {total_sum:.4e} / {4*math.pi:.4f}")
    print(f"  = {gamma_SI:.4e} Gy*m^2/(Bq*s)")

    # Step 4: Convert to uSv*m^2/(MBq*h)
    conv = 1e6 * 1e6 * 3600  # Sv->uSv * Bq->MBq * s->h
    gamma_final = gamma_SI * conv
    print(f"\n--- Step 4: Convert to uSv*m^2/(MBq*h) ---")
    print(f"  Gy->Sv: x1 (Q=1, photon)")
    print(f"  Sv->uSv: x1e6")
    print(f"  Bq->MBq: x1e6")
    print(f"  s->h: x3600")
    print(f"  Factor = 1e6 * 1e6 * 3600 = {conv:.2e}")
    print(f"  Gamma = {gamma_SI:.4e} * {conv:.2e} = {gamma_final:.4f} uSv*m^2/(MBq*h)")

    # Step 5: Convert to R*cm^2/(mCi*h) for comparison
    gamma_R = gamma_final / 0.023703
    print(f"\n--- Step 5: Reverse to R*cm^2/(mCi*h) ---")
    print(f"  = {gamma_final:.4f} / 0.023703 = {gamma_R:.3f} R*cm^2/(mCi*h)")

    # Comparison with Smith & Stabin
    if name == "Am-241":
        ss_R = 0.749
        ss_uSv = 0.01775
    else:
        ss_R = 0.858
        ss_uSv = 0.02034

    print(f"\n--- Comparison with Smith & Stabin ---")
    print(f"  NIST calc: {gamma_R:.3f} R*cm^2/(mCi*h) = {gamma_final:.4f} uSv*m^2/(MBq*h)")
    print(f"  S&S [R1]:  {ss_R:.3f} R*cm^2/(mCi*h) = {ss_uSv:.5f} uSv*m^2/(MBq*h)")
    print(f"  Ratio NIST/S&S: {gamma_R/ss_R:.4f}")
    print(f"  Difference: {(gamma_R/ss_R - 1)*100:+.1f}%")

    # Per-line summary table
    print(f"\n--- Per-line summary ---")
    print(f"{'Line':<12} {'E(keV)':>8} {'Y':>8} {'muen_air':>10} {'contrib':>12} {'frac':>8}")
    for E, Y, label, mu, E_J, mu_m2, c in contributions:
        frac = c / total_sum * 100
        print(f"{label:<12} {E:>8.2f} {Y:>8.4f} {mu:>10.4f} {c:>12.4e} {frac:>7.1f}%")
    print(f"{'TOTAL':<12} {'':>8} {'':>8} {'':>10} {total_sum:>12.4e} {'100.0%':>8}")
    print()
