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

# NIST Iron mu/rho (keV, cm^2/g) - above K-edge 7.112 keV
fe_data = [
    (10, 170.6), (15, 57.08), (20, 25.68), (30, 8.176),
    (40, 3.629), (50, 1.958), (60, 1.205), (80, 0.5952), (100, 0.3717)
]
rho_fe = 7.874

# Lead segments with edge handling
rho_pb = 11.35

def get_pb_mu(E):
    if E < 13.035:
        return log_interp(E, 10, 130.6, 13.035, 67.01)
    elif E < 15.2:
        return log_interp(E, 13.035, 162.1, 15.2, 107.8)
    elif E < 15.861:
        return log_interp(E, 15.2, 148.5, 15.861, 134.4)
    elif E <= 20:
        return log_interp(E, 15.861, 154.8, 20, 86.36)
    elif E <= 30:
        return log_interp(E, 20, 86.36, 30, 30.32)
    elif E <= 40:
        return log_interp(E, 30, 30.32, 40, 14.36)
    elif E <= 50:
        return log_interp(E, 40, 14.36, 50, 8.041)
    elif E <= 60:
        return log_interp(E, 50, 8.041, 60, 5.021)
    elif E <= 80:
        return log_interp(E, 60, 5.021, 80, 2.419)
    elif E < 88.005:
        return log_interp(E, 80, 2.419, 88.005, 1.910)
    elif E <= 100:
        return log_interp(E, 88.005, 7.683, 100, 5.549)
    else:
        return log_interp(E, 100, 5.549, 150, 2.014)

def get_fe_mu(E):
    # [Gemini 교차검증 수정] 원래 코드는 E < fe_data[0][0] 일 때 루프를 통과하여
    # fe_data[-1][1] (100 keV 값)을 반환하는 버그가 있었음.
    # 예: Fe-55 (5.89 keV) 입력 시 → μ/ρ=0.3717 (100 keV) 반환 → 차폐 1,916배 과소평가.
    # 수정: 배열 범위 이하 에너지에 대해 처음 두 점으로 log-log 외삽 적용.
    if E < fe_data[0][0]:
        e1, m1 = fe_data[0]
        e2, m2 = fe_data[1]
        return log_interp(E, e1, m1, e2, m2)
    for i in range(len(fe_data)-1):
        e1, m1 = fe_data[i]
        e2, m2 = fe_data[i+1]
        if e1 <= E <= e2:
            return log_interp(E, e1, m1, e2, m2)
    return fe_data[-1][1]

def get_air_muen(E):
    # [Gemini 교차검증 수정] 동일한 fallback 버그. 범위 이하 에너지 → 외삽 처리.
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

def calc_weights(lines):
    contribs = []
    for E, Y, name in lines:
        muen = get_air_muen(E)
        gi = Y * E * muen
        contribs.append((E, Y, name, muen, gi))
    total = sum(c[4] for c in contribs)
    return [(E, Y, name, muen, gi, gi/total) for E, Y, name, muen, gi in contribs], total

def transmission(x_cm, weights, mu_func, rho):
    T = 0
    for E, Y, name, muen, gi, wi in weights:
        mu = mu_func(E) * rho
        T += wi * math.exp(-mu * x_cm)
    return T

def find_thickness(weights, mu_func, rho, target_T):
    x_lo, x_hi = 0, 10
    while transmission(x_hi, weights, mu_func, rho) > target_T:
        x_hi *= 2
    for _ in range(100):
        x_mid = (x_lo + x_hi) / 2
        T = transmission(x_mid, weights, mu_func, rho)
        if abs(T - target_T) < 1e-10:
            break
        if T > target_T:
            x_lo = x_mid
        else:
            x_hi = x_mid
    return x_mid

print("=" * 80)
print("NUCLIDE COMPOSITE SPECTRUM HVL CALCULATION")
print("=" * 80)

for name, lines in [("Am-241", am241_lines), ("Cd-109", cd109_lines)]:
    print(f"\n{'='*80}")
    print(f"  {name}")
    print(f"{'='*80}")

    weights, total = calc_weights(lines)

    print(f"\n--- Emission line contributions ---")
    print(f"{'Line':<16} {'E(keV)':>8} {'Y':>8} {'muen_air':>10} {'Y*E*muen':>10} {'w_i':>8}")
    print("-" * 68)
    for E, Y, nm, muen, gi, wi in weights:
        print(f"{nm:<16} {E:>8.2f} {Y:>8.4f} {muen:>10.4f} {gi:>10.4f} {wi:>8.4f}")
    print(f"{'TOTAL':<16} {'':>8} {'':>8} {'':>10} {total:>10.4f} {'1.0000':>8}")

    # X-ray vs gamma fraction
    if name == "Am-241":
        xray_frac = sum(wi for E, Y, nm, muen, gi, wi in weights if E < 25)
        gamma_frac = sum(wi for E, Y, nm, muen, gi, wi in weights if E >= 25)
        print(f"\n  X-ray fraction (E<25keV): {xray_frac*100:.1f}%")
        print(f"  Gamma fraction (E>=25keV): {gamma_frac*100:.1f}%")
    elif name == "Cd-109":
        xray_frac = sum(wi for E, Y, nm, muen, gi, wi in weights if E < 30)
        gamma_frac = sum(wi for E, Y, nm, muen, gi, wi in weights if E >= 30)
        print(f"\n  X-ray fraction (E<30keV): {xray_frac*100:.1f}%")
        print(f"  Gamma fraction (E>=30keV): {gamma_frac*100:.1f}%")

    # Iron composite
    hvl_fe = find_thickness(weights, get_fe_mu, rho_fe, 0.5)
    qvl_fe = find_thickness(weights, get_fe_mu, rho_fe, 0.25)
    tvl_fe = find_thickness(weights, get_fe_mu, rho_fe, 0.1)

    print(f"\n--- Fe composite (rho=7.874) ---")
    print(f"  HVL = {hvl_fe*10:.4f} mm")
    print(f"  QVL = {qvl_fe*10:.4f} mm")
    print(f"  TVL = {tvl_fe*10:.4f} mm")

    # Lead composite
    hvl_pb = find_thickness(weights, get_pb_mu, rho_pb, 0.5)
    qvl_pb = find_thickness(weights, get_pb_mu, rho_pb, 0.25)
    tvl_pb = find_thickness(weights, get_pb_mu, rho_pb, 0.1)
    cvl_pb = find_thickness(weights, get_pb_mu, rho_pb, 0.01)
    mvl_pb = find_thickness(weights, get_pb_mu, rho_pb, 0.001)

    print(f"\n--- Pb composite (rho=11.35) ---")
    print(f"  HVL = {hvl_pb*10:.5f} mm")
    print(f"  QVL = {qvl_pb*10:.5f} mm")
    print(f"  TVL = {tvl_pb*10:.5f} mm")
    print(f"  CVL = {cvl_pb*10:.5f} mm")
    print(f"  MVL = {mvl_pb*10:.5f} mm")

    # Smith & Stabin comparison
    if name == "Am-241":
        ss = {"HVL": 0.00974, "QVL": 0.0235, "TVL": 0.106, "CVL": 0.528, "MVL": 0.948}
        nist = {"HVL": hvl_pb*10, "QVL": qvl_pb*10, "TVL": tvl_pb*10, "CVL": cvl_pb*10, "MVL": mvl_pb*10}
    else:
        ss = {"HVL": 0.0109, "QVL": 0.023, "TVL": 0.0353, "CVL": 0.0805, "MVL": 0.74}
        nist = {"HVL": hvl_pb*10, "QVL": qvl_pb*10, "TVL": tvl_pb*10, "CVL": cvl_pb*10, "MVL": mvl_pb*10}

    print(f"\n--- S&S vs NIST comparison (Pb) ---")
    print(f"  {'':>5} {'S&S(mm)':>10} {'NIST(mm)':>10} {'ratio':>8}")
    for k in ["HVL", "QVL", "TVL", "CVL", "MVL"]:
        ratio = nist[k] / ss[k] if ss[k] > 0 else float('inf')
        print(f"  {k:>5} {ss[k]:>10.5f} {nist[k]:>10.5f} {ratio:>8.3f}")

    # Transmission table
    print(f"\n--- Transmission vs thickness ---")
    for mat, mf, rho, thicknesses in [
        ("Fe", get_fe_mu, rho_fe, [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]),
        ("Pb", get_pb_mu, rho_pb, [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    ]:
        print(f"\n  {mat}:")
        for t_mm in thicknesses:
            T = transmission(t_mm/10, weights, mf, rho)
            if T > 0.001:
                print(f"    {t_mm:>7.3f} mm  T={T*100:>9.4f}%")
            else:
                print(f"    {t_mm:>7.3f} mm  T={T:.3e}")
