import math
import lab1
import matplotlib.pyplot as plt
import kursach as kurs
from weasyprint import HTML
import coeficents


#ПРЕДРАСЧЁТ
gC           =  0.86
sigmaBX      = 0.99
sigmaCC      = 0.955
sigma1       = 0.99
sigma2       = 0.955
phiC1        = 0.99
phiC2        = 0.97
effkGas      = 0.995
effkCompFull = 0.81
effkLpcFull = 0.87
effkFanFull  = 0.92
aCoeff       = 0.03
effkHptFull  = 0.91
effkLptFull  = 0.92
ksiTake      = 0.02
gAirBack     = 0.014
z = 0.5

cB   = 220.0
cBX  = cB - 45.0 #45
cKND = 180.0
cB2  = cB - 17.5 #17.5
cK = 155.0
cG   = 155.0 # 135
MHPT = 0.5 # 0.6
MLPT = 0.4 # 0.325
dB = 0.35
dHPC = 0.6
D_cp_H = 11

fBlade = (0.25 + 0.35) / 2
uB1 = (380 + 490) / 2
uK1 = (450 + 500) / 2
mu = (1.2 + 1.8) / 2
blade_knd_density = 4450 #BT-6
blade_kvd_density = 8063 #INCONEL 718
blade_tnd_density = 8080 #UDIMET710
blade_tvd_density = 8400 #ЖС6У
DcpDivHvyhHighPressure = (6.0 + 20.0) / 2
DcpDivHvyhOther = 3.0
F1           = (0.25 + 0.35) / 2
P0           = 101325.0
T0           = 288.0
ro0          = P0 / (287.0 * T0)

var_calc_B = 0 # 0-вт,1-ср,2-корпус
var_calc_rzd = 0
var_calc_knd = 1
var_calc_vkvd = 1
var_calc_K = 1
var_calc_G = 0
var_calc_tnd = 0
#ПРЕДРАСЧЁТ
expected_P = 264447 * 1.07
m = 5.3 # 
TGfull = 1701.5428
TKfull = 838.2633856500527
PiK_full = 34.7929
variant = lab1.calc_proto(coeficents.coef,TGfull, m,PiK_full)
TKfull = variant.T_k
alpha = variant.alpha
Xopt = variant.x_opt
k_pred = variant.k
k_gas_pred = variant.k_gas
P_spec_pred = variant.p_spec
c_spec_pred = variant.c_spec
Lfree = variant.l_free_energy
print("TG* =", TGfull,"\n"
    "m =", m,"\n"
    "pik* =", PiK_full,"\n"
    "TK* =", TKfull,"\n"
    "alpha =", alpha,"\n"
    "Xopt =", Xopt,"\n"
    "k =", k_pred,"\n"
    "k' =", k_gas_pred,"\n"
    "P_spec =", P_spec_pred,"\n"
    "c_spec =", c_spec_pred,"\n"
    "Lfree =", Lfree,"\n"
    )


def calc_density(P, T, R):
    return P / (T * R)
def calculate_LB2(q_T, v_or6,x, L_CB, eta_THA_star, m):
    return ((1 + q_T - v_or6) * x * L_CB * eta_THA_star) / m
T0 = 288.0
P0 = 101325.0


def main():
#region Сечения
#Сечение ВХ-ВХ
    print("BX-ВX:")
    TBXfull = T0
    PBXfull = P0
    TBX = TBXfull - (cBX * cBX / 2 / 1004.5)
    PBX = PBXfull * math.pow(TBXfull/TBX, (1.4/(1-1.4)))
    roBX = PBX / (287 * TBX)
    print("\tTBX* =", TBXfull,
        "\n\tTBX =", TBX,
        "\n\tPBX* =", PBXfull,
        "\n\tPBX =", PBX,
        "\n\troBX =", roBX)

#Сечение В-В
    print("В-В:")
    PBfull = PBXfull * sigmaBX
    TBfull = TBXfull
    cP = 1004.5
    TB = TBfull - (cB * cB / 2 /cP)
    PB = PBfull * math.pow(TBfull/TB, (1.4/(1-1.4)))
    roB = calc_density(PB, TB, 287.0)
    print("\tTB* =", TBfull,
        "\n\tTB =", TB,
        "\n\tPB* =", PBfull,
        "\n\tPB =", PB,
        "\n\troB =", roB)
    
#Доп параметры для сечений
    print("Params:")
    compressor = kurs.calculate_compressor_temperature(
                288, PiK_full, coeficents.coef["effk_comp_full"]
            )
    combustion_props = kurs.calculate_combustion_properties(
                0.86, 0.14, TKfull, TGfull, coeficents.coef["effk_gas"]
            )
    q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
    print("\tq_T = ",q_T)
    Hu = combustion_props["Hu"]
    v_take = coeficents.coef["ksi_take"] - coeficents.coef["g_air_back"]
    LB2 = calculate_LB2(q_T,v_take,Xopt,Lfree,coeficents.coef["effk_lpt_full"],6.7)
    print("\tLB2 = ",LB2)
    #TB2full = kurs.calculate_fan_temperature(TBfull,LB2,coeficents.coef["effk_fan_full"])[0]
    PiBstar = kurs.calculate_fan_temperature(TBfull,LB2,coeficents.coef["effk_fan_full"])[1]
    #print("\tTB2full = ",TB2full)
    print("\tPiBfull = ",PiBstar)
    LB2 = kurs.calc_LB2(TBfull,PiBstar,coeficents.coef["effk_fan_full"])["LB2"]
    TB2full = kurs.calc_LB2(TBfull,PiBstar,coeficents.coef["effk_fan_full"])["TB2full"]
    print("\tLB2 = ", LB2,"\n\tTB2full = ",TB2full)
    #Принял итоговое пи вентилятора =1,7 т.к при пи=1.95 полный пи всего КНД получался 1,42, что бред
    
    LK = kurs.calc_LK(TBfull, TKfull, PiK_full, effkCompFull)["LK"]
    print("\tLK =",LK)
    
    LKND = (z * (LK * effkCompFull + m * LB2 * effkFanFull) - m * LB2 * effkFanFull) / effkLpcFull
    print("\tLKND =", LKND)
    
    PiKND_full = kurs.calc_PiKND(TBfull, LKND, effkLpcFull)["PiKND_full"]
    print("\tPiKND* =",PiKND_full)
    TKNDfull = kurs.calc_PiKND(TBfull, LKND, effkLpcFull)["TKNDfull"]
    print("\tTKND* =",TKNDfull) 
    
    PiKVD_full = PiK_full / PiKND_full
    print("\tPiKVD* =",PiKVD_full)
    
    LKVD_adiabatic = kurs.calc_LKVD(TKNDfull, TKfull, PiKVD_full)["LKVD_adiabatic"]
    print("\tLKVD adiabatic =",LKVD_adiabatic)
    
    effkHpcFull = LKVD_adiabatic / (LK - LKND)
    print("\teffk HPC* =", effkHpcFull)

#ПЕРВЫЙ КОНТУР

#Cечение КНД-КВД
    print("KНД-KВД:")
    cpLPC = kurs.get_cp_real(TKNDfull)["cp"]
    kLPC = kurs.get_cp_real(TKNDfull)["k"]
    print("\n\tcpLPC* =",cpLPC,"\n\tkLPC* =",kLPC)
    PKNDfull = PBfull * PiKND_full
    print("\tPKND* =",PKNDfull)
    TKND = TKNDfull - cKND *cKND / 2 / cpLPC
    print("\tTKND =",TKND)
    PKND = PKNDfull * math.pow(TKNDfull / TKND, kLPC / (1 - kLPC))
    print("\tPKND =",PKND)
    roKND = calc_density(PKND, TKND, 287.0)
    print("\troKND =",roKND)

#Сечение К-К
    print("K-K:")
    cpHPC = kurs.get_cp_real(TKfull)["cp"]
    kHPC = kurs.get_cp_real(TKfull)["k"]
    print("\tcpHPC* =",cpHPC,"\n\tkHPC* =",kHPC)
    PKfull = PKNDfull * PiKVD_full
    print("\tPK* =",PKfull)
    LKVD = LKVD_adiabatic / effkHpcFull
    print("\tLKVD =",LKVD)
    TK = TKfull - cK *cK / 2 / cpHPC
    print("\tTK =",TK)
    PK = PKfull * math.pow(TKfull / TK, kHPC / (1 - kHPC))
    print("\tPK =",PK)
    roK = calc_density(PK, TK, 287.0)
    print("\troK =",roK)
#Сечение Г-Г
    print("Г-Г:")
    cpGas = kurs.get_cpgas(0.86, 0.14, TGfull, alpha)["cp_mix"]
    kGas = kurs.get_cpgas(0.86, 0.14, TGfull, alpha)["k"]
    RGas = kurs.get_cpgas(0.86, 0.14, TGfull, alpha)["R_mix"]
    print("\tcpGas* =",cpGas,"\n\tkGas* =",kGas,"\n\tRGas* =",RGas)
    PGfull = PKfull * sigmaCC
    print("\tPG* =",PGfull)
    TG = TGfull - cG *cG / 2 / cpGas
    print("\tTG =",TG)
    PG = PGfull * math.pow(TGfull / TG, kGas / (1 - kGas))
    print("\tPG =",PG)
    roG = calc_density(PG, TG, 287.0)
    print("\troG =",roG)
    #уточнение охлаждения, коррекция температуры Тг*
    TGcorr = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["TGcorr"]
    cpMix = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
    kMix = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["k"]
    RMix = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    print("\tTGorr = ",TGcorr,"\n\tcpMix = ", cpMix,"\n\tkMix = ",kMix,"\n\tRMix = ", RMix)
    
#Сечение ТВД-ТВД
    print("ТВД-ТВД:")
    LTVD = LKVD / (1 + q_T - v_take)
    print("\tLTVD =",LTVD)
    TTVDfull = kurs.calc_LTurb(TGcorr, alpha, LTVD, effkHptFull, ksiTake, q_T, gAirBack)["TTurb"]
    print("\tTTVD* =",TTVDfull)
    PiTVD_full = kurs.calc_LTurb(TGcorr, alpha, LTVD, effkHptFull, ksiTake, q_T, gAirBack)["PiTurb"]
    print("\tPiTVD* =",PiTVD_full)
    PTVDfull = PGfull / PiTVD_full
    print("\tPTVD* =", PTVDfull)
    kMixTVDfull = kurs.get_cpmix(TTVDfull, alpha, ksiTake, q_T, gAirBack)["k"]
    RMixTVDfull = kurs.get_cpmix(TTVDfull, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    print("\tkMixTVDfull = ",kMixTVDfull,"\n\tRMixTVDfull = ", RMixTVDfull)
    TTVD = TTVDfull / (1 + (kMixTVDfull - 1) * MHPT * MHPT / 2)
    print("\tTTVD =", TTVD)
    PTVD = PTVDfull * math.pow(TTVDfull / TTVD, kMixTVDfull / (1 - kMixTVDfull))
    print("\tPTVD =", PTVD)
    kMixTVD = kurs.get_cpmix(TTVD, alpha, ksiTake, q_T, gAirBack)["k"]
    RMixTVD = kurs.get_cpmix(TTVD, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    roTVD = calc_density(PTVD, TTVD, RMixTVD)
    print("\troTVD =",roTVD)
    print("\tkMixTVD = ",kMixTVD,"\n\tRMixTVD = ", RMixTVD)
    cTVD = MHPT * math.sqrt(kMixTVD * RMixTVD * TTVD)
    print("\tcTVD =", cTVD)
    
#Сечение ТНД-ТНД
    print("ТНД-ТНД:")
    LTND = (LB2 * m + LKND) / (1 + q_T - v_take)
    print("\tLTND =", LTND)
    TTNDfull = kurs.calc_LTurb(TTVDfull, alpha, LTND, effkLptFull, ksiTake, q_T, gAirBack)["TTurb"]
    print("\tTTND* =",TTNDfull)
    PiTND_full = kurs.calc_LTurb(TTVDfull, alpha, LTND, effkLptFull, ksiTake, q_T, gAirBack)["PiTurb"]
    print("\tPiTND* =",PiTND_full)
    PTNDfull = PTVDfull / PiTND_full
    print("\tPTND* =", PTNDfull)
    kMixTNDfull = kurs.get_cpmix(TTNDfull, alpha, ksiTake, q_T, gAirBack)["k"]
    RMixTNDfull = kurs.get_cpmix(TTNDfull, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    print("\tkMixTVDfull = ",kMixTVDfull,"\n\tRMixTVDfull = ", RMixTVDfull)
    TTND = TTNDfull / (1 + (kMixTNDfull - 1) * MLPT * MLPT / 2)
    print("\tTTND =", TTND)
    PTND = PTNDfull * math.pow(TTNDfull / TTND, kMixTNDfull / (1 - kMixTNDfull))
    print("\tPTND =", PTND)
    kMixTND = kurs.get_cpmix(TTND, alpha, ksiTake, q_T, gAirBack)["k"]
    RMixTND = kurs.get_cpmix(TTND, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    roTND = calc_density(PTND, TTND, RMixTND)
    print("\troTND =",roTND)
    print("\tkMixTVD = ",kMixTVD,"\n\tRMixTVD = ", RMixTVD)
    cTND = MLPT * math.sqrt(kMixTND * RMixTND * TTND)
    print("\tcTND =", cTND)
    
#Сечение С1-С1
    print("С1-С1:")
    TC1full = TTNDfull
    print("\tTC1* =", TC1full)
    kMixC1 = kurs.get_cpmix(TC1full, alpha, ksiTake, q_T, gAirBack)["k"]
    RMixC1 = kurs.get_cpmix(TC1full, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    kMixC1ave = kurs.get_cpmix_ave(T0, TC1full, alpha, ksiTake, q_T, gAirBack)["k"]
    cpMixC1ave = kurs.get_cpmix_ave(T0, TC1full, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
    RMixC1ave = kurs.get_cpmix_ave(T0, TC1full, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    print("\tkMixC1 = ",kMixC1)
    cc1t = math.sqrt(2 * kMixC1ave / (kMixC1ave + 1) * RMixC1ave * TC1full)
    print("\tcc1t =", cc1t)
    cc1 = phiC1 * cc1t
    print("\tcc1 =", cc1)
    sigmaC1 = kurs.GDF_pressure(kMixC1ave, 1, 1)["GDF_P"] / kurs.GDF_pressure(kMixC1ave, 1, phiC1)["GDF_P"]
    GDF1 = kurs.GDF_pressure(kMixC1ave, 1, phiC1)["GDF_P"]
    print("\tGDF1 =", GDF1)
    print("\tsigmaC1 =", sigmaC1)
    PiC1_full = PTNDfull * sigma1 / P0
    print("\tPiC1* =", PiC1_full)
    PiCr = math.pow((kMixC1 + 1) / 2, kMixC1 / (kMixC1 - 1))
    print("\tPiCr* =", PiCr)
    cc1_docr = kurs.calc_cc1_docr(TC1full, T0, PiC1_full, phiC1, alpha, ksiTake, q_T, gAirBack)["cc1_docr"]
    print("\tcc1 =", cc1_docr)
    PC1full = PTNDfull * sigmaC1
    print("\tPC1* =", PC1full)
    PC1 = P0
    print("\tPC1 =", PC1)

    print("\tkMixC1ave = ",kMixC1ave,"\n\tcpMixC1ave = ", cpMixC1ave,"\n\tRMixC1ave = ", RMixC1ave)
    TC1 = TC1full - cc1_docr * cc1_docr / 2 / cpMixC1ave
    print("\tTC1 =", TC1)
    roC1 = calc_density(PC1, TC1, RMixC1ave)
    print("\troC1 =",roC1)
#ВТОРОЙ КОНТУР

#Сечение В2-В2
    print("В2-В2:")
    PB2full = PBfull * PiBstar
    print("\tPB2* =", PB2full)
    cpB2 = kurs.get_cp_real(TB2full)["cp"]
    kB2 = kurs.get_cp_real(TB2full)["k"]
    TB2 = TB2full - cB2 * cB2 / 2 / cpB2
    print("\tTB2 =", TB2)
    PB2 = PB2full * math.pow(TB2full / TB2, kB2 / (1 - kB2))
    print("\tPB2 =", PB2)
    roB2 = calc_density(PB2, TB2, 287.0)
    print("\troB2 =", roB2)
    
#Сечение С2-С2
    print("С2-С2:")
    TC2full = TB2full
    print("\tTC2* =", TC2full)
    PC2full = PB2full * sigma2
    print("\tPC2* =", PC2full)
    kC2 = kurs.get_cp_real(TC2full)["k"]
    print("\tkC2 = ",kC2)
    sigmaC2 = kurs.GDF_pressure(kC2, 1, 1)["GDF_P"] / kurs.GDF_pressure(kC2, 1, phiC2)["GDF_P"]
    print("\tsigmaC2 =", sigmaC2) #дать комментарий почему так
    GDF2a = kurs.GDF_pressure(kC2, 1, 1)["GDF_P"]
    GDF2b = kurs.GDF_pressure(kC2, 1, phiC2)["GDF_P"]
    print("\tGDF2a = ",GDF2a,"\n\tGDF2b = ", GDF2b)
    PiC2_full = PC2full / P0 / sigmaC2
    print("\tPiC2* =", PiC2_full)
    Pi2Cr = math.pow((kC2 + 1) / 2, kC2 / (kC2 - 1))
    print("\tPi2Cr* =", Pi2Cr)
    cpC2 = kurs.get_cp_air(TC2full, T0)["cp"]
    kC2 = kurs.get_cp_air(TC2full, T0)["k"]
    cc2_docr = kurs.calc_cc2_docr(TC2full, T0, PiC2_full, phiC2)["cc2_docr"]
    print("\tcc2 =", cc2_docr)
    PC2 = P0
    print("\tPC2 =", PC2)
    TC2 = TC2full - cc2_docr * cc2_docr / 2 / cpC2
    print("\tTC2 =", TC2)
    roC2 = calc_density(PC2, TC2, 287.0)
    print("\troC2 =", roC2)
    
#Основные параметры двигателя
    print("Основные параметры двигателя:")
    P_spec = cc1_docr * (1 + q_T - v_take) / (m + 1) + cc2_docr * m / (m + 1)
    print("\tP spec =", P_spec)
    Gair = expected_P / P_spec
    print("\tGair =", Gair)
    Gair1 = Gair / (1 + m)
    print("\tGair1 =", Gair1)
    Gair2 = Gair - Gair1
    print("\tGair2 =", Gair2)
    Ggas = (1 + q_T - v_take) * Gair1
    print("\tGgas =", Ggas)
    C_spec = (3600 * q_T * (1 - ksiTake)) / ((1 + m) * P_spec)
    print("\tC_spec =", C_spec)
    EngEff = ((1 + q_T - v_take) * (cc1_docr * cc1_docr / 2) + (m * cc2_docr * cc2_docr / 2)) / (Hu * q_T * effkGas)
    print("\tEngine eff =", EngEff)
    NB2 = LB2 * Gair2
    print("\tNB2 =", NB2)
    NKND = LKND * Gair1
    print("\tNKND =", NKND)
    NKVD = LKVD * Gair1
    print("\tNKVD =", NKVD)
    NTVD = LTVD * Ggas
    print("\tNTVD =", NTVD)
    NTND = LTND * Ggas
    print("\tNTND =", NTND)
    P_diff = (P_spec_pred - P_spec) / P_spec_pred * 100
    print("\tP diff = ",P_diff)
    Ndiff = (NKND + NB2 - NTND)
    print("\tN diff = ",Ndiff)
#endregion Cечения

#region Геометрия
    #Геометрия
    # Исходные параметры (заполни нужными значениями)
    F_BX = Gair / (cBX * roBX)
    # Диаметр круглого входного сечения:
    D_BX = math.sqrt(4 * F_BX / math.pi)
    
        # Исходные параметры
# Вентилятор
    Dcp_B = 0.0
    Dvt_B = 0.0
    d_vt_B = dB
    F_B = Gair / (cB * roB)  
    if var_calc_B == 0:
        Dvt_B = math.sqrt(4 * F_B / (math.pi * (1 /(d_vt_B * d_vt_B) - 1)))
        D_B = Dvt_B / d_vt_B
        Dcp_B = (D_B + Dvt_B) / 2
    elif var_calc_B == 1:
    # относительный диаметр втулки вентилятора
        Dcp_B = math.sqrt((F_B * (1 + d_vt_B)) / (math.pi * (1 - d_vt_B)))
        # Диаметр корпуса и втулки вентилятора
        D_B = (2 * Dcp_B) / (1 + d_vt_B)
        Dvt_B = D_B * d_vt_B
    elif var_calc_B == 2:
        D_B = math.sqrt(4 * F_B / (math.pi * (1 - d_vt_B * d_vt_B)))
        Dvt_B = D_B * d_vt_B
        Dcp_B = (D_B + Dvt_B) / 2
    h_B = (D_B - Dvt_B) / 2
    
#Разделитель
    Dcp_rzd = 0.0
    Dvt_rzd = 0.0
    D_rzd = 0.0
    F_rzd = Gair / (cB2 * roB2)
    var_calc_rzd = var_calc_B
    if var_calc_rzd == 0:
        Dvt_rzd = Dvt_B
        D_rzd = math.sqrt(Dvt_rzd ** 2 + (4 * F_rzd / math.pi))
        Dcp_rzd = (D_rzd + Dvt_rzd) / 2
    elif var_calc_rzd == 1:
        Dcp_rzd = Dcp_B
        Dvt_rzd = Dcp_rzd - (F_rzd / (math.pi * Dcp_rzd))
        D_rzd = 2 * Dcp_rzd - Dvt_rzd
    elif var_calc_rzd == 2:
        D_rzd = D_B
        Dvt_rzd = math.sqrt(D_rzd ** 2 - (4 * F_rzd / math.pi))
        Dcp_rzd = (D_rzd + Dvt_rzd) / 2
    # Высота лопатки СА вентилятора (или мнимая высота)
    h_SA_rzd = (D_rzd - Dvt_rzd) / 2
    #внутренний
    #
    # Площадь сечения канала второго контура
    F_B2 = Gair2 / (cB2 * roB2)
    # Диаметр разделителя
    D_vrzd = math.sqrt(D_rzd**2 - (4 * F_B2 / math.pi))
    # Высота канала второго контура или высота лопатки НА вентилятора
    h_B2 = (D_rzd - D_vrzd) / 2
    # Толщина разделителя
    b_rzd = 0.1 * h_SA_rzd  # ... или 0.1 * h_SA_red
    # Внешний диаметр канала первого контура
    F1 = Gair1 / (cB2 * roB2)
    D_vknd = D_vrzd - 2 * b_rzd
    # Втулочный диаметр канала первого контура
    D_vt_vknd = math.sqrt(D_vknd ** 2 - (4 * F1 / math.pi))
    # Средний диаметр канала первого контура
    D_cp_vknd = (D_vknd + D_vt_vknd) / 2
    # Высота лопатки на входе в канал первого контура
    h_vknd = (D_vknd - D_vt_vknd) / 2
    # Диаметр корпуса и втулки на выходе из КНД
#КНД
    D_knd = 0.0
    Dvt_knd = 0.0
    Dcp_knd = 0.0
    F_knd = Gair1 / (cKND * roKND)    
    if var_calc_knd == 0:
        Dvt_knd = D_vt_vknd
        D_knd = math.sqrt(Dvt_knd ** 2 + (4 * F_knd / math.pi))
        Dcp_knd = (D_knd + Dvt_knd) / 2
    elif var_calc_knd == 1:
        Dcp_knd = D_cp_vknd
        Dvt_knd = Dcp_knd - (F_knd / math.pi / Dcp_knd)
        D_knd = (2 * Dcp_knd - Dvt_knd)
    elif var_calc_knd == 2:
        D_knd = D_vknd
        Dvt_knd = math.sqrt(D_knd ** 2 - (4 * F_knd / math.pi))
        Dcp_knd = (D_knd + Dvt_knd) / 2
    # Высота лопатки последнего СА КНД
    h_knd = (D_knd - Dvt_knd) / 2
    # Проверка допустимости геометрии КНД
    # 1. Проверка минимальной длины лопатки
    min_blade_length = 0.012  # 12 мм в метрах
    is_min_length_ok = h_knd >= min_blade_length
    # 2. Проверка относительного диаметра втулки
    relative_hub_diameter = Dvt_knd / D_knd
    max_relative_hub_diameter = 0.92  # максимально допустимое значение
    min_relative_hub_diameter = 0.82  # минимальное рекомендуемое значение
    is_hub_diameter_ok = (relative_hub_diameter <= 0.82)
    

    # Общая проверка геометрии КНД
    is_knd_geometry_valid = is_min_length_ok and is_hub_diameter_ok
#вКВД-вКВД

    D_vkvd = 0.0
    Dvt_vkvd = 0.0
    Dcp_vkvd = 0.0
    F_vkvd = F_knd #Так как пи не равны
    d_vt_kvd = dHPC #!!!!!!!!!!!!!НЕУВЕРЕН
    if var_calc_vkvd == 0:
        Dvt_vkvd = math.sqrt((4 * F_vkvd) / (math.pi * (1 / (d_vt_kvd ** 2) - 1)))
        D_vkvd = Dvt_vkvd / d_vt_kvd
        Dcp_vkvd = (D_vkvd + Dvt_vkvd) / 2
    elif var_calc_vkvd == 1:
        Dcp_vkvd = math.sqrt((F_vkvd * (1 + d_vt_kvd)) / (math.pi * (1 - d_vt_kvd)))
        D_vkvd = (2 * Dcp_vkvd) / (1 + d_vt_kvd)
        Dvt_vkvd = D_vkvd * d_vt_kvd
    elif var_calc_vkvd == 2:
        D_vkvd = math.sqrt((4 * F_vkvd) / (math.pi * (1 - d_vt_kvd * d_vt_kvd)))
        Dvt_vkvd = D_vkvd * d_vt_kvd
        Dcp_vkvd = (D_vkvd + Dvt_vkvd) / 2
    
    h_vkvd = (D_vkvd - Dvt_vkvd) / 2
#K-K
    # Площадь сечения на выходе из КВД

    F_K = Gair1 / (cK * roK) #!!!!!!!!!! почему cK а не cK вот roK норм
    D_K = 0.0
    Dvt_K = 0.0
    Dcp_K = 0.0
    if var_calc_K == 0:
        Dvt_K = Dvt_vkvd
        D_K = math.sqrt(Dvt_K**2 + (4 * F_K / math.pi))
        Dcp_K = (D_K + Dvt_K) / 2
    elif var_calc_K == 1:
        Dcp_K = Dcp_vkvd
        Dvt_K = Dcp_K - (F_K / math.pi / Dcp_K )
        D_K = (2 * Dcp_K) - Dvt_K
    elif var_calc_K == 2:
        D_K = D_vkvd
        Dvt_K = math.sqrt(D_K**2 - (4 * F_K / math.pi))
        Dcp_K = (D_K + Dvt_K) / 2
    h_K = (D_K - Dvt_K) / 2
    #Проверка
    compressor_blade_check = h_K >= 0.015
    rel_vt_vkvd_check = (Dvt_vkvd / D_vkvd) >= 0.5
    rel_vt_K_check = 0.92 >= (Dvt_K / D_K) >= 0.87
    rel_blade_check = 6.25 >= (h_vknd / h_K) >= 2
    compressor_check = (compressor_blade_check
                        and rel_vt_vkvd_check
                        and rel_vt_K_check
                        and rel_blade_check)
##ТВД-ТВД
    #Dср твд / h твд = 13
    F_tvd = Gair1 * (1 + q_T - v_take) / (cTVD * roTVD)
    h_tvd = math.sqrt((4 * F_tvd) / (math.pi*(math.pow((D_cp_H + 1),2) - math.pow((D_cp_H - 1),2))))
    D_tvd = h_tvd * (D_cp_H + 1)
    Dvt_tvd = h_tvd * (D_cp_H - 1)
    Dcp_tvd = D_cp_H * h_tvd
    # Площадь проходного сечения Г-Г (без учета возврата воздуха)
#Г-Г
    D_g = 0.0
    Dcp_g = 0.0
    Dvt_g = 0.0
    F_g = (Gair1 * (1 + q_T - ksiTake)) / (cG * roG)
    if var_calc_G == 0:
        Dvt_g = Dvt_tvd
        D_g = math.sqrt(Dvt_g ** 2 + (4 * F_g / math.pi))
        Dcp_g = (D_g + Dvt_g) / 2
    elif var_calc_G == 1:
        Dcp_g = Dcp_tvd
        Dvt_g = Dcp_g - (F_g / math.pi / Dcp_g)
        D_g = 2 * Dcp_g - Dvt_g
    elif var_calc_G == 2:
        D_g = D_tvd
        Dvt_g = math.sqrt(D_g ** 2 - (4 * F_g / math.pi))
        Dcp_g = (D_g + Dvt_g) / 2
    h_g = (D_g - Dvt_g) / 2
    #Проверка ТВД
    tvd_blade_check = 5.9 >= (h_tvd / h_g) >= 1.1
#ТНД-ТНД
    F_tnd = Gair1 * (1 + q_T - v_take) / (cTND * roTND)
    Dcp_tnd = 0.0
    Dvt_tnd = 0.0
    D_tnd = 0.0
    if var_calc_tnd == 0:
        Dvt_tnd = Dvt_tvd
        D_tnd = math.sqrt(Dvt_tnd ** 2 + (4 * F_tnd / math.pi))
        Dcp_tnd = (D_tnd + Dvt_tnd) / 2
    elif var_calc_tnd == 1:
        Dcp_tnd = Dcp_tvd
        Dvt_tnd = Dcp_tnd - (F_tnd  / math.pi / Dcp_tnd)
        D_tnd = 2 * Dcp_tnd - Dvt_tnd
    elif var_calc_tnd == 2:
        D_tnd = D_tvd
        Dvt_tnd= math.sqrt(D_tnd ** 2 - (4 * F_tnd / math.pi))
        Dcp_tnd = (D_tnd + Dvt_tnd) / 2 
    h_tnd = (D_tnd - Dvt_tnd) / 2
    #Проверка ТНД
    tnd_blade_check = 5.9 >= (h_tnd / h_tvd) >= 1.1
    rel_blade_tnd = 7.5 >= Dcp_tnd / h_tnd >= 2.7
    turbine_check = tvd_blade_check and tnd_blade_check and rel_blade_check  
# C1-C1
    F_c1 = (Gair1 * (1 + q_T - ksiTake)) / (cc1 * roC1)
    D_c1 = math.sqrt(4 * F_c1 / math.pi)
    #вС2 надо сначала посчитать все остальное
#C2-C2
    D_vc2 = D_B
    Dvt_vc2 = D_rzd
    if D_vrzd < D_g or D_vrzd < D_tvd or D_vrzd < D_tnd: #Затуп, но в нашем случае это и не надо
        D_max = max(D_vknd,D_knd,D_vkvd,D_K,D_g,D_tvd,D_tnd)
        Dvt_vc2 = D_max
        
    F_c2 = Gair2 / (cc2_docr * roC2)
    Dvt_c2 = Dvt_vc2
    D_c2 = math.sqrt((4 * F_c2) / math.pi + Dvt_c2)
    Dcp_c2 = (D_c2 + Dvt_c2) / 2
#endregion Геометрия
#region Частоты ротторов и кол-во ступеней ТНД
    ucp_B = kurs.calc_u_cp(uB1,Dcp_B,D_B)
    ucp_vKND = kurs.calc_u_cp(uB1,D_cp_vknd,D_B)
    ucp_KND = kurs.calc_u_cp(uB1,Dcp_knd,D_B)
    ucp_TND = kurs.calc_u_cp(uB1,Dcp_tnd,D_B)
    ucp_TVD = (0.533*math.pow((1+m),0.536) + 0.6) * ucp_TND
    ucp_vKVD = kurs.calc_u_cp(ucp_TVD,Dcp_vkvd,Dcp_tvd)
    ucp_KVD = kurs.calc_u_cp(ucp_TVD,Dcp_K,Dcp_tvd)
    ucp_G = kurs.calc_u_cp(ucp_TVD,Dcp_g,Dcp_tvd)
    ucp_TVD = kurs.calc_u_cp(ucp_TVD,Dcp_tvd,Dcp_tvd)
    uk1_check = (ucp_TVD * D_vkvd / Dcp_tvd) < 500
    if uk1_check:
        print(f"Окружная скорость на среднем диаметре входа в компрессор в норме [{ucp_vKVD}]")
    else:
        print(f"[ОШИБКА] Окружная скорость на среднем диаметре входа в компрессор не в норме [{ucp_vKVD}]")
    sigma_B_KND = 600 * 10 ** 6
    sigma_B_KVD = 724 * 10 ** 6
    sigma_B_TVD = 380 * 10 ** 6
    sigma_B_TND = 475 * 10 ** 6
    sigma_p_KND = kurs.calc_sigma_p(blade_knd_density, ucp_KND, Dvt_knd / D_knd, fBlade)
    sigma_p_KVD = kurs.calc_sigma_p(blade_kvd_density, ucp_KVD, Dvt_K / D_K, fBlade)
    sigma_p_TVD = kurs.calc_sigma_p(blade_tvd_density, ucp_TVD, Dvt_tvd / D_tvd, fBlade)
    sigma_p_TND = kurs.calc_sigma_p(blade_tnd_density, ucp_TND, Dvt_tnd / D_tnd, fBlade)
    kurs.check_sigma(sigma_p_KND, sigma_B_KND)
    kurs.check_sigma(sigma_p_KVD, sigma_B_KVD)
    kurs.check_sigma(sigma_p_TVD, sigma_B_TVD)   
    kurs.check_sigma(sigma_p_TND, sigma_B_TND)
    n_VD = ucp_TVD / (math.pi * Dcp_tvd) * 60 #min
    n_ND = ucp_TND / (math.pi * Dcp_tnd) * 60 #min
    turbine_res = kurs.calc_stages(LTVD,effkHptFull,Dcp_tvd,Dcp_g,ucp_TVD,uK1,D_vkvd)
    stages_TVD = turbine_res["z"]
    y = turbine_res["y"]
    mu_tvd =  turbine_res["mu"]
    print(f"Количество ступеней ТВД: {stages_TVD}, с запасом {y*100}%, коэф нагрузки {mu_tvd}")
    turbine_res_ND = kurs.calc_stages(LTND,effkLptFull,Dcp_tnd,Dcp_tvd,ucp_TND,uB1,D_vknd)
    stages_TND = turbine_res_ND["z"]
    y_ND = turbine_res_ND["y"]
    mu_tnd =  turbine_res_ND["mu"]
    print(f"Количество ступеней ТНД: {stages_TND}, с запасом {y_ND*100}%, коэф нагрузки {mu_tnd}")

#endregion Частоты ротторов и кол-во ступеней ТНД

#region График сечений
# Построение второго графика на правой оси
    temp = [1 * 2,2 * 2,3 * 2,5 * 2,10 * 2,3 * 2,4 * 2,7 * 2,3 * 2]
    x = formArray(temp)
    labels = ['Н',"ВХ","В","КНД","К","Г","ТВД","ТНД","С1"]
    
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8*2.5, 6*2.5))
    first_axes = axes[0]
    first_axes.set_xticks(x)
    first_axes.set_xticklabels(labels, rotation=0, ha='center')
    y1 = [T0,TBXfull,TBfull,TKNDfull,TKfull,TGfull,TTVDfull,TTNDfull,TC1full]
    y2 = [T0,TBX,TB,TKND,TK,TG,TTVD,TTND,TC1]
    first_axes.plot(x, y1, label='Полная', color='red')
    first_axes.legend(loc='best')
    first_axes.set_ylabel('Температуры')
    first_axes.tick_params(axis='y', labelcolor='b')
    first_axes.set_ylim(200,1800)
    right_ax = first_axes.twinx()
    right_ax.plot(x, y2, label='Статика', color='red', linestyle='-.')
    lines1, labels1 = first_axes.get_legend_handles_labels()
    lines2, labels2 = right_ax.get_legend_handles_labels()

    first_axes.legend(lines1 + lines2, labels1 + labels2, loc='best')

    right_ax.set_ylim(200,1800)
    for x_i in x:
        plt.axvline(x=x_i,color='black', linestyle='--')
        
    second_axes = axes[1]
    y1 = [P0,PBXfull,PBfull,PKNDfull,PKfull,PGfull,PTVDfull,PTNDfull,PC1full]
    y2 = [P0,PBX,PB,PKND,PK,PG,PTVD,PTND,PC1]
    second_axes.set_xticks(x)
    second_axes.set_xticklabels(labels, rotation=0, ha='center')
    second_axes.plot(x, y1,label = 'Полная',color ='blue')
    second_axes.set_ylabel('Давления')
    second_axes.tick_params(axis='y', labelcolor='b')
    second_axes.set_ylim(20*1000,3750*1000)
    right_ax_sec = second_axes.twinx()
    right_ax_sec.plot(x, y2, label = 'Статика',color = 'blue',linestyle='-.')
    right_ax_sec.set_ylim(20*1000,3750*1000)
    
    lines1, labels1 = second_axes.get_legend_handles_labels()
    lines2, labels2 = right_ax_sec.get_legend_handles_labels()

    second_axes.legend(lines1 + lines2, labels1 + labels2, loc='best')
    for x_i in x:
        plt.axvline(x=x_i, color='black', linestyle='--')
        
    thir_axes = axes[2]
    y2 = [0,cBX,cB,cKND,cK,cG,cTVD,cTND,cc1]
    thir_axes.set_xticks(x)
    thir_axes.set_xticklabels(labels, rotation=0, ha='center')
    thir_axes.plot(x, y2, 'black')
    thir_axes.set_ylabel('Скорости')
    thir_axes.tick_params(axis='y', labelcolor='b')
    second_axes.legend(lines1 + lines2, labels1 + labels2, loc='best')
    for x_i in x:
        thir_axes.axvline(x=x_i, color='black', linestyle='--')
    plt.tight_layout() # Автоматически корректирует расположение

    # 4. Отображение графика
    plt.savefig('../diagrams/sections.png', dpi=300, bbox_inches='tight')
#endregion График сечений
#region График Геометрии
    figNew, axesNew = plt.subplots(nrows=1, ncols=1, figsize=(26, 13))
    first_axes_new = axesNew
    temp = [0,1,2,3,5,6,10,12,13,15,17]
    x = temp
    labels = ['Вх','В',"Рзд","КНДв","КНДвых","КВДвх","КВДвых","Г","ТВД","ТНД","С1"]
    first_axes_new.set_xticks(x)
    first_axes_new.set_xticklabels(labels, rotation=0, ha='center')
    Dvt = [0,Dvt_B / 2 ,Dvt_rzd / 2,D_vt_vknd / 2,Dvt_knd / 2,Dvt_vkvd / 2,Dvt_K / 2,Dvt_g / 2,Dvt_tvd / 2,Dvt_tnd / 2,0]
    Dcp = [D_BX / 4,Dcp_B / 2,Dcp_rzd / 2,D_cp_vknd / 2,Dcp_knd / 2,Dcp_vkvd / 2,Dcp_K / 2,Dcp_g / 2,Dcp_tvd / 2,Dcp_tnd / 2,D_c1 / 4]
    D = [D_BX / 2,D_B / 2,D_rzd / 2,D_vknd / 2,D_knd / 2,D_vkvd / 2,D_K / 2,D_g / 2,D_tvd / 2,D_tnd / 2,D_c1 / 2]
    y_soplo_second = [D_BX / 2,D_B / 2,D_B / 2,D_B / 2,D_B / 2,D_vc2 / 2,D_c2 / 2,D_c2 / 2,D_c2 / 2,D_c2 / 2,D_c2 / 4]
    first_axes_new.plot(x, Dvt, label='Втулка', color='black')
    first_axes_new.plot(x, Dcp, label='Средняя', color='orange', linestyle='-.')
    first_axes_new.plot(x, D, label='Корпус', color='black')
    first_axes_new.plot(x, y_soplo_second, label='С2', color='blue')
    first_axes_new.legend(loc='best')
    for x_i in x:
        plt.axvline(x=x_i,color='black', linestyle='--')
    plt.savefig("../diagrams/comp.png")
#endregion График Геометрии
#region Вывод
    generate_engine_report_html({
                "Коэффициенты": [
            ("g_c", gC),
            ("sigma_intake", sigmaBX),
            ("sigma_cc", sigmaCC),
            ("sigma_1", sigma1),
            ("sigma_2", sigma2),
            ("phi_c1", phiC1),
            ("phi_c2", phiC2),
            ("effk_gas", effkGas),
            ("effk_comp_full", effkCompFull),
            ("effk_fan_full", effkFanFull),
            ("a", aCoeff),
            ("effk_hpt_full", effkHptFull),
            ("effk_lpt_full", effkLptFull),
            ("ksi_take", ksiTake),
            ("g_air_back", gAirBack),
            ("z", z)
        ],

        "Скорости и диаметры": [
            ("cB", cB),
            ("cBX", cBX),        # cB - 45
            ("cKND", cKND),
            ("cB2", cB2),        # cB - 16
            ("cK", cK),
            ("cG", cG),
            ("MHPT", MHPT),
            ("MLPT", MLPT),
            ("dB", dB),          # (0.30 + 0.65)/2
            ("dHPC", dHPC)         # (0.50 + 0.65)/2
        ],

        "Прочие параметры": [
            ("DcpDivHvyhHighPressure",D_cp_H),  # (6.0 + 20.0)/2
            ("DcpDivHvyhOther", 3.0),
            ("F1", F1),                        # (0.25 + 0.35)/2
            ("P0", P0),
            ("T0", T0),
            ("ro0", ro0)                       # P0 / (287.0 * T0)
        ],

        "Предрасчёт": [
            ("expected_P", expected_P),          # 264447 * 1.07
            ("m", m),
            ("TGfull", TGfull),
            ("TKfull", TKfull),
            ("PiK_full", PiK_full)
        ],
        "Params": [
            ("q_T", q_T),
            ("LB2", LB2),
            ("TB2*", TB2full),
            ("PiB*", PiBstar),
            ("LK", LK),
            ("LKND", LKND),
            ("PiKND*", PiKND_full),
            ("TKND*", TKNDfull),
            ("PiKVD*", PiKVD_full),
            ("LKVD_adiabatic", LKVD_adiabatic),
            ("effk_HPC*", effkHpcFull),
        ],
        "BX-BX": [
            ("TBX*", TBXfull),
            ("TBX", TBX),
            ("PBX*", PBXfull),
            ("PBX", PBX),
            ("roBX", roBX),
        ],
        "B-B": [
            ("TB*", TBfull),
            ("TB", TB),
            ("PB*", PBfull),
            ("PB", PB),
            ("roB", roB),
        ],
        "KНД-KВД": [
            ("cpLPC*", cpLPC),
            ("kLPC*", kLPC),
            ("PKND*", PKNDfull),
            ("TKND", TKND),
            ("PKND", PKND),
            ("roKND", roKND),
        ],
        "К-К" : [
            ("cpHPC*", cpHPC),
            ("kHPC*", kHPC),
            ("PK*", PKfull),
            ("LKVD", LKVD),
            ("TK", TK),
            ("PK", PK),
            ("roK", roK),
        ],
        "Г-Г": [
            ("cpGas*", cpGas),
            ("kGas*", kGas),
            ("RGas*", RGas),
            ("PG*", PGfull),
            ("TG", TG),
            ("PG", PG),
            ("roG", roG),
            ("tTGorr", TGcorr),
            ("tcpMix", cpMix),
            ("kMix", kMix),
            ("RMix", RMix)
        ],
        "TВД-ТВД": [
            ("LTVD", LTVD),
            ("TTVD*", TTVDfull),
            ("PiTVD*", PiTVD_full),
            ("PTVD*", PTVDfull),
            ("kMixTVDfull", kMixTVDfull),
            ("RMixTVDfull", RMixTVDfull),
            ("TTVD", TTVD),
            ("PTVD", PTVD),
            ("kMixTVD", kMixTVD),
            ("RMixTVD", RMixTVD),
            ("roTVD", roTVD),
            ("cTVD", cTVD)
        ],
        "ТНД-ТНД": [
            ("LTND", LTND),
            ("TTND*", TTNDfull),
            ("PiTND*", PiTND_full),
            ("PTND*", PTNDfull),
            ("kMixTNDfull", kMixTNDfull),
            ("RMixTNDfull", RMixTNDfull),
            ("TTND", TTND),
            ("PTND", PTND),
            ("kMixTND", kMixTND),
            ("RMixTND", RMixTND),
            ("roTND", roTND),
            ("cTND", cTND)
        ],

        "С1-С1": [
            ("TC1*", TC1full),
            ("kMixC1", kMixC1),
            ("RMixC1", RMixC1),
            ("cc1t", cc1t),
            ("cc1", cc1),
            ("sigmaC1", sigmaC1),
            ("GDF1", GDF1),
            ("PiC1*", PiC1_full),
            ("PiCr*", PiCr),
            ("cc1_docr", cc1_docr),
            ("PC1*", PC1full),
            ("PC1", PC1),
            ("kMixC1ave", kMixC1ave),
            ("cpMixC1ave", cpMixC1ave),
            ("RMixC1ave", RMixC1ave),
            ("TC1", TC1),
            ("roC1", roC1)
        ],

        "В2-В2": [
            ("PB2*", PB2full),
            ("cpB2", cpB2),
            ("kB2", kB2),
            ("TB2", TB2),
            ("PB2", PB2),
            ("roB2", roB2)
        ],

        "С2-С2": [
            ("TC2*", TC2full),
            ("PC2*", PC2full),
            ("kC2", kC2),
            ("sigmaC2", sigmaC2),
            ("GDF2a", GDF2a),
            ("GDF2b", GDF2b),
            ("PiC2*", PiC2_full),
            ("Pi2Cr*", Pi2Cr),
            ("cpC2", cpC2),
            ("TC2", TC2),
            ("PC2", PC2),
            ("roC2", roC2)
        ],

        "Основные параметры двигателя": [
            ("P_spec", P_spec),
            ("Gair", Gair),
            ("Gair1", Gair1),
            ("Gair2", Gair2),
            ("Ggas", Ggas),
            ("C_spec", C_spec),
            ("EngEff", EngEff),
            ("NB2", NB2),
            ("NKND", NKND),
            ("NKVD", NKVD),
            ("NTVD", NTVD),
            ("NTND", NTND),
            ("P_diff", P_diff),
            ("Ndiff", Ndiff)
        ],
        "Окружные скорости": [
        ("Вход в вентилятор(В)", ucp_B),
        ("Вход в КНД(вКНД)", ucp_vKND),
        ("Выход КНД(КНД)", ucp_KND),
        ("Вход в КВД(вКВД)", ucp_vKVD),
        ("Выход КВД(КВД)", ucp_KVD),
        ("Вход в ТВД(Г)", ucp_G),
        ("Выход ТВД(ТВД)", ucp_TVD),
        ("Выход ТНД(ТНД)", ucp_TND),
        ("Частота вращения ротора ВД", n_VD),
        ("Частота вращения ротора НД", n_ND),
    ],

    "прочность": [
        (" σб KND", sigma_B_KND),
        (" σб KVD", sigma_B_KVD),
        (" σб TVD", sigma_B_TVD),
        (" σб TND", sigma_B_TND),
        ("σp KND", sigma_p_KND),
        ("σp KVD", sigma_p_KVD),
        ("σp TVD", sigma_p_TVD),
        ("σp TND", sigma_p_TND),
    ],

    "Ступени ТВД": [
        ("Кол-во ступеней", stages_TVD),
        ("параметр Парсонса", y),
    ],

    "Ступени ТНД": [
        ("Кол-во ступеней", stages_TND),
        ("параметр Парсонса", y_ND),
    ],
            
},"../calcs/sections_report.html", "../calcs/sections_report.pdf")
    def printGeometry():
        print("=" * 50)
        print("РАСЧЕТ ГЕОМЕТРИИ ПРОТОЧНОЙ ЧАСТИ ГТД")
        print("=" * 50)

        print("\n" + "=" * 30)
        print("СЕЧЕНИЕ ВХОДА (BX)")
        print("=" * 30)
        print(f"Площадь входного сечения F_BX: {F_BX:.4f} м²")
        print(f"Диаметр входного сечения D_BX: {D_BX:.4f} м")

        print("\n" + "=" * 30)
        print("СЕЧЕНИЕ ВЕНТИЛЯТОРА (B)")
        print("=" * 30)
        print(f"Площадь тракта вентилятора F_B: {F_B:.4f} м²")
        print(f"Относительный диаметр втулки d_vt_B: {d_vt_B:.3f}")
        print(f"Средний диаметр Dcp_B: {Dcp_B:.4f} м")
        print(f"Диаметр корпуса D_B: {D_B:.4f} м")
        print(f"Диаметр втулки Dvt_B: {Dvt_B:.4f} м")
        print(f"Высота первой лопатки h_B: {h_B:.4f} м")

        print("\n" + "=" * 30)
        print("СЕЧЕНИЕ РАЗДЕЛИТЕЛЯ (P3J)")
        print("=" * 30)
        print(f"Площадь перед разделителем F_rzd: {F_rzd:.4f} м²")
        print(f"Средний диаметр Dcp_rzd: {Dcp_rzd:.4f} м")
        print(f"Диаметр корпуса D_rzd: {D_rzd:.4f} м")
        print(f"Диаметр втулки Dvt_rzd: {Dvt_rzd:.4f} м")
        print(f"Высота лопатки СА h_SA_rzd: {h_SA_rzd:.4f} м")
        print(f"Толщина разделителя b_rzd: {b_rzd:.4f} м")

        print("\n" + "=" * 30)
        print("ВТОРОЙ КОНТУР (B2)")
        print("=" * 30)
        print(f"Площадь канала F_B2: {F_B2:.4f} м²")
        print(f"Диаметр разделителя D_rzd: {D_rzd:.4f} м")
        print(f"Высота канала h_B2: {h_B2:.4f} м")

        print("\n" + "=" * 30)
        print("КНД - ВХОД (vKND)")
        print("=" * 30)
        print(f"Высота лопатки h_vknd: {h_vknd:.4f} м")
        print(f"Внешний диаметр D_vknd: {D_vknd:.4f} м")
        print(f"Втулочный диаметр D_vt_vknd: {D_vt_vknd:.4f} м")
        print(f"Средний диаметр D_cp_vknd: {D_cp_vknd:.4f} м")

        print("\n" + "=" * 30)
        print("КНД - ВЫХОД")
        print("=" * 30)
        print(f"Площадь сечения F_knd: {F_knd:.4f} м²")
        print(f"Средний диаметр Dcp_knd: {Dcp_knd:.4f} м")
        print(f"Диаметр корпуса D_knd: {D_knd:.4f} м")
        print(f"Диаметр втулки Dvt_knd: {Dvt_knd:.4f} м")
        print(f"Высота лопатки h_knd: {h_knd:.4f} м")
        print(f"Относительный диаметр втулки: {relative_hub_diameter:.3f}")
        print(f"Геометрия КНД корректна: {'ДА' if is_knd_geometry_valid else 'НЕТ'}")
        
        print("\n" + "=" * 30)
        print("КВД - ВХОД (vKVD)")
        print("=" * 30)
        print(f"Относительный диаметр втулки d_vt_kvd: {d_vt_kvd:.3f}")
        print(f"Внешний диаметр D_vkvd: {D_vkvd:.4f} м")
        print(f"Диаметр втулки Dvt_vkvd: {Dvt_vkvd:.4f} м")
        print(f"Средний диаметр Dcp_vkvd: {Dcp_vkvd:.4f} м")
        print(f"Высота лопатки h_vkvd: {h_vkvd:.4f} м")

        print("\n" + "=" * 30)
        print("КВД - ВЫХОД (K)")
        print("=" * 30)
        print(f"Площадь сечения F_K: {F_K:.4f} м²")
        print(f"Внешний диаметр D_K: {D_K:.4f} м")
        print(f"Диаметр втулки Dvt_K: {Dvt_K:.4f} м")
        print(f"Средний диаметр Dcp_K: {Dcp_K:.4f} м")
        print(f"Высота лопатки h_K: {h_K:.4f} м")

        print("\n" + "=" * 30)
        print("ТВД")
        print("=" * 30)
        print(f"Площадь сечения F_tvd: {F_tvd:.4f} м²")
        print(f"Отношение D_cp_H: {D_cp_H}")
        print(f"Высота лопатки h_tvd: {h_tvd:.4f} м")
        print(f"Внешний диаметр D_tvd: {D_tvd:.4f} м")
        print(f"Диаметр втулки Dvt_tvd: {Dvt_tvd:.4f} м")
        print(f"Средний диаметр Dcp_tvd: {Dcp_tvd:.4f} м")

        print("\n" + "=" * 30)
        print("СЕЧЕНИЕ Г-Г")
        print("=" * 30)
        print(f"Площадь сечения F_g: {F_g:.4f} м²")
        print(f"Внешний диаметр D_g: {D_g:.4f} м")
        print(f"Диаметр втулки Dvt_g: {Dvt_g:.4f} м")
        print(f"Средний диаметр Dcp_g: {Dcp_g:.4f} м")
        print(f"Высота лопатки h_g: {h_g:.4f} м")

        print("\n" + "=" * 30)
        print("ТНД")
        print("=" * 30)
        print(f"Площадь сечения F_tnd: {F_tnd:.4f} м²")
        print(f"Средний диаметр Dcp_tnd: {Dcp_tnd:.4f} м")
        print(f"Диаметр втулки Dvt_tnd: {Dvt_tnd:.4f} м")
        print(f"Внешний диаметр D_tnd: {D_tnd:.4f} м")
        print(f"Высота лопатки h_tnd: {h_tnd:.4f} м")
        print(f"Геометрия ТНД корректна: {'ДА' if turbine_check else 'НЕТ'}")

        print("\n" + "=" * 30)
        print("СОПЛО ПЕРВОГО КОНТУРА (C1)")
        print("=" * 30)
        print(f"Площадь сечения F_c1: {F_c1:.4f} м²")
        print(f"Диаметр сопла D_c1: {D_c1:.4f} м")

        print("\n" + "=" * 30)
        print("ВХОД В СОПЛО ВОРОГО КОНТУРА (C2)")
        print("=" * 30)
        print(f"Площадь сечения F_c2: {F_c2:.4f} м²")
        print(f"Диаметр сопла D_vc2: {D_vc2:.4f} м")
        print(f"Внутрений диаметр сопла D_vc2: {Dvt_vc2:.4f} м")
        
        print("\n" + "=" * 30)
        print("СОПЛО ВОРОГО КОНТУРА (C2)")
        print("=" * 30)
        print(f"Площадь сечения F_c2: {F_c2:.4f} м²")
        print(f"Диаметр сопла D_c2: {D_c2:.4f} м")
        print(f"Внутрений диаметр сопла D_c2: {Dvt_c2:.4f} м")
        
        print("\n" + "=" * 50)
        print("ПРОВЕРКИ ГЕОМЕТРИИ")
        print("=" * 50)
        print(f"КНД - мин. высота лопатки ({min_blade_length*1000:.0f} мм): {'СОБЛЮДЕНО' if is_min_length_ok else 'НЕ СОБЛЮДЕНО'}")
        print(f"КНД - отн. диаметр втулки ({min_relative_hub_diameter}-{max_relative_hub_diameter}): {'СОБЛЮДЕНО' if is_hub_diameter_ok else 'НЕ СОБЛЮДЕНО' }")
        if compressor_check:
            print(f"ГЕОМЕТРИЯ КОМПРЕССОРА: СОБЛЮДЕНА")
        else:
            print(f"КОМПРЕССОР ЧЕК:\nВысота лопатки последнего СА:{"СОБЛЮДЕНО" if compressor_blade_check else "НЕ СОБЛЮДЕНО"}")
            print(f"Относительный диаметр втулки в сечении К:{"СОБЛЮДЕНО" if rel_vt_K_check else "НЕ СОБЛЮДЕНО - "}{Dvt_vkvd / D_vkvd}")
            print(f"Относительный диаметр втулки в сечении вКВД:{"СОБЛЮДЕНО" if rel_vt_vkvd_check else "НЕ СОБЛЮДЕНО"}")
            print(f"отношение длин лопаток:{"СОБЛЮДЕНО" if rel_blade_check else "НЕ СОБЛЮДЕНО"}")
        if turbine_check:
            print(f"ГЕОМЕТРИЯ ТУРБИНЫ: СОБЛЮДЕНА")
        else:
            print(f"ТУРБИНА ЧЕК:\nОтношение длин лопаток ТВД:{"СОБЛЮДЕНО" if tvd_blade_check else "НЕ СОБЛЮДЕНО"}")
            print(f"Отношение длин лопаток ТНД:{"СОБЛЮДЕНО" if rel_blade_tnd else "НЕ СОБЛЮДЕНО"}")
            print(f"Приведенная длина лопатки в сечении ТНД:{"СОБЛЮДЕНО" if rel_blade_check else "НЕ СОБЛЮДЕНО"}")
            
    def save_geometry_html(filename="geometry_report.html", pdf_filename="geometry_report.pdf"):
        def make_section(title, rows):
            html = f"<h2>{title}</h2>\n"
            html += "<table>\n<tr><th>Параметр</th><th>Значение</th></tr>\n"
            for name, value in rows:
                html += f"<tr><td>{name}</td><td>{value}</td></tr>\n"
            html += "</table>\n"
            return html

        # ---- CSS ----
        css = """
        <style>
            body {
                font-family: Arial, sans-serif;
                padding: 15px;
                background: #fafafa;
            }
            h2 {
                margin-top: 30px;
                color: #333;
                border-bottom: 2px solid #4a90e2;
                padding-bottom: 5px;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin-top: 10px;
                margin-bottom: 30px;
            }
            table, th, td {
                border: 1px solid #aaa;
            }
            th {
                background: #e8f0ff;
                font-weight: bold;
            }
            td, th {
                padding: 8px 12px;
                text-align: left;
            }
            tr:nth-child(even) td {
                background: #f6f6f6;
            }
            .graph-block {{ margin-top: 10px; text-align: center; }}
        </style>
        """

        # ---- HTML сборка ----
        html = "<html><head><meta charset='utf-8'>" + css + "</head><body>"
        html += "<h1>Отчёт по геометрии двигателя</h1>"

        html += make_section("СЕЧЕНИЕ ВХОДА (BX)", [
            ("Площадь входного сечения F_BX", f"{F_BX:.4f} м²"),
            ("Диаметр входного сечения D_BX", f"{D_BX:.4f} м"),
        ])

        html += make_section("СЕЧЕНИЕ ВЕНТИЛЯТОРА (B)", [
            ("Площадь тракта вентилятора F_B", f"{F_B:.4f} м²"),
            ("Относительный диаметр втулки d_vt_B", f"{d_vt_B:.3f}"),
            ("Средний диаметр Dcp_B", f"{Dcp_B:.4f} м"),
            ("Диаметр корпуса D_B", f"{D_B:.4f} м"),
            ("Диаметр втулки Dvt_B", f"{Dvt_B:.4f} м"),
            ("Высота первой лопатки h_B", f"{h_B:.4f} м"),
        ])

        html += make_section("СЕЧЕНИЕ РАЗДЕЛИТЕЛЯ (РЗД)", [
            ("Площадь перед разделителем F_rzd", f"{F_rzd:.4f} м²"),
            ("Средний диаметр Dcp_rzd", f"{Dcp_rzd:.4f} м"),
            ("Диаметр корпуса D_rzd", f"{D_rzd:.4f} м"),
            ("Диаметр втулки Dvt_rzd", f"{Dvt_rzd:.4f} м"),
            ("Высота лопатки СА h_SA_rzd", f"{h_SA_rzd:.4f} м"),
            ("Толщина разделителя b_rzd", f"{b_rzd:.4f} м"),
        ])

        html += make_section("ВТОРОЙ КОНТУР (B2)", [
            ("Площадь канала F_B2", f"{F_B2:.4f} м²"),
            ("Диаметр разделителя D_rzd", f"{D_rzd:.4f} м"),
            ("Высота канала h_B2", f"{h_B2:.4f} м"),
        ])

        html += make_section("КНД – ВХОД (vKND)", [
            ("Высота лопатки h_vknd", f"{h_vknd:.4f} м"),
            ("Внешний диаметр D_vknd", f"{D_vknd:.4f} м"),
            ("Втулочный диаметр D_vt_vknd", f"{D_vt_vknd:.4f} м"),
            ("Средний диаметр D_cp_vknd", f"{D_cp_vknd:.4f} м"),
        ])

        html += make_section("КНД – ВЫХОД", [
            ("Площадь сечения F_knd", f"{F_knd:.4f} м²"),
            ("Средний диаметр Dcp_knd", f"{Dcp_knd:.4f} м"),
            ("Диаметр корпуса D_knd", f"{D_knd:.4f} м"),
            ("Диаметр втулки Dvt_knd", f"{Dvt_knd:.4f} м"),
            ("Высота лопатки h_knd", f"{h_knd:.4f} м"),
            ("Относительный диаметр втулки", f"{relative_hub_diameter:.3f}"),
            ("Геометрия КНД корректна", "ДА" if is_knd_geometry_valid else "НЕТ"),
        ])

        html += make_section("КВД – ВХОД (vKVD)", [
            ("Относительный диаметр втулки d_vt_kvd", f"{d_vt_kvd:.3f}"),
            ("Внешний диаметр D_vkvd", f"{D_vkvd:.4f} м"),
            ("Диаметр втулки Dvt_vkvd", f"{Dvt_vkvd:.4f} м"),
            ("Средний диаметр Dcp_vkvd", f"{Dcp_vkvd:.4f} м"),
            ("Высота лопатки h_vkvd", f"{h_vkvd:.4f} м"),
        ])

        html += make_section("КВД – ВЫХОД", [
            ("Площадь сечения F_K", f"{F_K:.4f} м²"),
            ("Внешний диаметр D_K", f"{D_K:.4f} м"),
            ("Диаметр втулки Dvt_K", f"{Dvt_K:.4f} м"),
            ("Средний диаметр Dcp_K", f"{Dcp_K:.4f} м"),
            ("Высота лопатки h_K", f"{h_K:.4f} м"),
        ])
        html += make_section("СЕЧЕНИЕ Г-Г", [
            ("Площадь сечения F_g", f"{F_g:.4f} м²"),
            ("Внешний диаметр D_g", f"{D_g:.4f} м"),
            ("Диаметр втулки Dvt_g", f"{Dvt_g:.4f} м"),
            ("Средний диаметр Dcp_g", f"{Dcp_g:.4f} м"),
            ("Высота лопатки h_g", f"{h_g:.4f} м"),
        ])
        html += make_section("ТВД", [
            ("Площадь сечения F_tvd", f"{F_tvd:.4f} м²"),
            ("Отношение D_cp_H", f"{D_cp_H}"),
            ("Высота лопатки h_tvd", f"{h_tvd:.4f} м"),
            ("Внешний диаметр D_tvd", f"{D_tvd:.4f} м"),
            ("Диаметр втулки Dvt_tvd", f"{Dvt_tvd:.4f} м"),
            ("Средний диаметр Dcp_tvd", f"{Dcp_tvd:.4f} м"),
        ])
        html += make_section("ТНД", [
            ("Площадь сечения F_tnd", f"{F_tnd:.4f} м²"),
            ("Средний диаметр Dcp_tnd", f"{Dcp_tnd:.4f} м"),
            ("Диаметр втулки Dvt_tnd", f"{Dvt_tnd:.4f} м"),
            ("Внешний диаметр D_tnd", f"{D_tnd:.4f} м"),
            ("Высота лопатки h_tnd", f"{h_tnd:.4f} м"),
            ("Геометрия ТНД корректна", "ДА" if turbine_check else "НЕТ"),
        ])
        html += make_section("СОПЛО ПЕРВОГО КОНТУРА", [
            ("Площадь сечения F_tnd", f"{F_c1:.4f} м²"),
            ("Диаметр сопла D_c1", f"{D_c1:.4f} м"),
        ])
        html += make_section("СОПЛО ВТОРОГО КОНТУРА", [
            ("Площадь сечения F_c2", f"{F_c2:.4f} м²"),
            ("Диаметр сопла D_c2", f"{D_c2:.4f} м"),
            ("Внутрений диаметр сопла Dvt_c2", f"{Dvt_c2:.4f} м"),
        ])
        html += make_section("ПРОВЕРКИ ГЕОМЕТРИИ", [
            ("КНД — мин. высота лопатки", "СОБЛЮДЕНО" if is_min_length_ok else "НЕ СОБЛЮДЕНО"),
            ("КНД — относительный диаметр втулки", "СОБЛЮДЕНО" if is_hub_diameter_ok else "НЕ СОБЛЮДЕНО"),
            ("Геометрия компрессора", "СОБЛЮДЕНА" if compressor_check else "НЕТ"),
            ("Геометрия турбины", "СОБЛЮДЕНА" if turbine_check else "НЕТ"),
        ])
        html += "</body></html>"

        # ---- Сохранить HTML ----
        with open(filename, "w", encoding="utf-8") as f:
            f.write(html)

        print(f"[OK] HTML файл сохранён: {filename}")

        # ---- Автоматическая генерация PDF ----
        try:
            HTML(filename).write_pdf(pdf_filename)
            print(f"[OK] PDF файл сохранён: {pdf_filename}")
        except Exception as e:
            print("[WARN] PDF не был создан:", e)
    printGeometry()
    save_geometry_html("../calcs/geometry_report.html", "../calcs/geometry.pdf")
#endregion Вывод


def formArray(arr):
    res = []
    temp = 0
    for value in arr:
        res.append(temp + value)
        temp += value
    return res
def generate_engine_report_html(data: dict, file_name_html,file_name_pdf) -> str:
    """
    Принимает словарь вида:
    {
        "BX-БX": [
            ("TBX*", TBXfull),
            ("TBX", TBX),
            ("PBX*", PBXfull),
            ("PBX", PBX),
            ("roBX", roBX),
        ],
        "B-B": [
            ("TB*", TBfull),
            ("TB", TB),
            ("PB*", PBfull),
            ("PB", PB),
            ("roB", roB),
        ],
        ...
    }
    И возвращает HTML со всеми таблицами.
    """

    # ----- CSS для красивых таблиц -----
    css = """
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background: #fafafa;
        }
        h2 {
            margin-top: 40px;
            color: #333;
            border-left: 6px solid #007acc;
            padding-left: 10px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
            margin-bottom: 30px;
            background: white;
        }
        th, td {
            border: 1px solid #bbb;
            padding: 8px 12px;
            text-align: left;
        }
        th {
            background: #e0e0e0;
        }
        tr:nth-child(even) {
            background: #f5f5f5;
        }
    </style>
    """

    html = f"<html><head><meta charset='utf-8'> {css} </head><body>"
    html += f"<h1>Расчёт параметров двигателя</h1>"

    # ---- Генерация таблиц ----
    for section, rows in data.items():
        html += f"<h2>{section}</h2>\n"
        html += "<table>\n"
        html +="<tr><th>Параметр</th><th>Значение</th></tr>\n"
        for name, value in rows:
            html += f"<tr><td>{name}</td><td>{value}</td></tr>\n"
        html += "</table>\n"

    html += "</body></html>"
    with open(file_name_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"[OK] HTML файл сохранён: {file_name_html}")

        # ---- Автоматическая генерация PDF ----
    try:
        HTML(file_name_html).write_pdf(file_name_pdf)
        print(f"[OK] PDF файл сохранён: {file_name_pdf}")
    except Exception as e:
        print("[WARN] PDF не был создан:", e)

# Вызов функции main
if __name__ == "__main__":
    main()