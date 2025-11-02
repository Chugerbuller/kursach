import math
import lab1
import kursach as kurs
coef = {
        "g_c": 0.86,
        "sigma_intake": 0.99,
        "sigma_cc": 0.955,
        "sigma_1": 0.99,
        "phi_c1": 0.98,
        "phi_c2": 0.98,
        "effk_gas": 0.995,
        "effk_comp_full": 0.81,
        "effk_fan_full": 0.87,
        "a": 0.03,
        "effk_hpt_full": 0.88,
        "effk_lpt_full": 0.91,
        "ksi_take": 0.155,
        "g_air_back": 0.124,
    }

cB   = 220.0
cBX  = cB - 45.0
cLPC = 180.0
cB2  = cB - 17.5
cHPC = 155.0
cG   = 135.0
MHPT = 0.6
MLPT = 0.325
dB = (0.30 + 0.65) / 2
dHPC = (0.50 + 0.65) / 2


DcpDivHvyhHighPressure = (6.0 + 20.0) / 2
DcpDivHvyhOther = 3.0
F1           = (0.25 + 0.35) / 2
P0           = 101325.0
T0           = 288.0
ro0          = P0 / (287.0 * T0)
gC           = 0.86
sigmaBX      = 0.99
sigmaCC      = 0.955
sigma1       = 0.99
phiC1        = 0.98
phiC2        = 0.98
effkGas      = 0.995
effkCompFull = 0.81
effkFanFull  = 0.87
aCoeff       = 0.03
effkHptFull  = 0.88
effkLptFull  = 0.91
ksiTake      = 0.155
gAirBack     = 0.124

def calc_density(P, T, R):
    return P / (T * R)
def calculate_LB2(q_T, v_or6,x, L_CB, eta_THA_star, m):
    return ((1 + q_T - v_or6) * x * L_CB * eta_THA_star) / m
T0 = 208.0
P0 = 101325.0
cBX = 175.0
cB = 220.0
sigmaBX = 0.99
m = 5.1
T = 1721.8024
PI_full = 27.6000
variant = lab1.calc_proto(coef,T, m,PI_full)
def main():
    TBXfull = T0
    PBXfull = P0
    TBX = TBXfull - (cBX * cBX / 2 / 1004.5)
    PBX = PBXfull * math.pow(TBXfull/TBX, (1.4/(1-1.4)))
    PoBX = PBX / (287 * TBX)
    print(f"{TBX:.2f}|{PBX:.2f}|{PoBX:.2f}")
    
    PBfull = PBXfull * sigmaBX
    TBfull = TBXfull
    cP = 1004.5
    TB = TBfull - (cB * cB / cP)
    PB = PBfull * math.pow(TBfull/TB, (1.4/(1-1.4)))
    roB = calc_density(PB, TB, 287.0)
    compressor = kurs.calculate_compressor_temperature(
                288, PI_full, coef["effk_comp_full"]
            )
    combustion_props = kurs.calculate_combustion_properties(
                0.86, 0.14, compressor['TK'], T, coef["effk_gas"]
            )
    q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
    v_take = coef["ksi_take"] - coef["g_air_back"]
    LB2 = calculate_LB2(q_T,v_take,variant.x_opt,variant.l_free_energy,coef["effk_lpt_full"],m)
    TB2full = kurs.calculate_fan_temperature(TBfull,LB2,coef["effk_fan_full"])[0]
    print(TB2full)
    LB2 = kurs.calc_LB2(TBfull,TB2full,1.95,coef["effk_fan_full"])["LB2"]
    TB2full = kurs.calc_LB2(TBfull,TB2full,1.95,coef["effk_fan_full"])["TB2full"]
    print(LB2,TB2full)

# Вызов функции main
if __name__ == "__main__":
    main()