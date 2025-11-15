import math
import lab1
import kursach as kurs
import matplotlib.pyplot as plt
import numpy as np
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
sigma2 = 0.99 #назначил больше диапазона
phiC1        = 0.98
phiC2        = 0.98
effkGas      = 0.995
effkCompFull = 0.81
effkLpcFull = 0.86
effkFanFull  = 0.87
aCoeff       = 0.03
effkHptFull  = 0.88
effkLptFull  = 0.91
ksiTake      = 0.155
gAirBack     = 0.124
z = 0.5
alpha = 2.5186

def calc_density(P, T, R):
    return P / (T * R)
def calculate_LB2(q_T, v_or6,x, L_CB, eta_THA_star, m):
    return ((1 + q_T - v_or6) * x * L_CB * eta_THA_star) / m
T0 = 288.0
P0 = 101325.0
cBX = 175.0
cB = 220.0
sigmaBX = 0.99
m = 4.9200
TGfull = 1721.8024
TKfull = 816.3413215872308
PiK_full = 27.6000
Xopt = 0.8022
Lfree = 578020.3908
variant = lab1.calc_proto(coef,TGfull, m,PiK_full)
def main():
    
#Сечение ВХ-ВХ
    TBXfull = T0
    PBXfull = P0
    TBX = TBXfull - (cBX * cBX / 2 / 1004.5)
    PBX = PBXfull * math.pow(TBXfull/TBX, (1.4/(1-1.4)))
    roBX = PBX / (287 * TBX)
    print("TBX* =", TBXfull,
          "TBX =", TBX,
          "PBX* =", PBXfull,
          "PBX =", PBX,
          "roBX =", roBX,
        
    )
    #print(f"{TBX:.2f}|{PBX:.2f}|{PoBX:.2f}")

#Сечение В-В
    PBfull = PBXfull * sigmaBX
    TBfull = TBXfull
    cP = 1004.5
    TB = TBfull - (cB * cB / cP)
    PB = PBfull * math.pow(TBfull/TB, (1.4/(1-1.4)))
    roB = calc_density(PB, TB, 287.0)
    print("TBX* =", TBfull,
          "TBX =", TB,
          "PBX* =", PBfull,
          "PBX =", PB,
          "roBX =", roB,
        
    )
    
#Доп параметры для сечений
    compressor = kurs.calculate_compressor_temperature(
                288, PiK_full, coef["effk_comp_full"]
            )
    combustion_props = kurs.calculate_combustion_properties(
                0.86, 0.14, compressor['TK'], TGfull, coef["effk_gas"]
            )
    q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
    Hu = combustion_props["Hu"]
    v_take = coef["ksi_take"] - coef["g_air_back"]
    LB2 = calculate_LB2(q_T,v_take,Xopt,Lfree,coef["effk_lpt_full"],6.7)
    print(LB2)
    TB2full = kurs.calculate_fan_temperature(TBfull,LB2,coef["effk_fan_full"])[0]
    print(TB2full)
    LB2 = kurs.calc_LB2(TBfull,1.7,coef["effk_fan_full"])["LB2"]
    TB2full = kurs.calc_LB2(TBfull,1.7,coef["effk_fan_full"])["TB2full"]
    print(LB2,TB2full)
    #Принял итоговое пи вентилятора =1,7 т.к при пи=1.95 полный пи всего КНД получался 1,42, что бред
    
    LK = kurs.calc_LK(TBfull, TKfull, PiK_full, effkCompFull)["LK"]
    print("LK =",LK)
    
    LKND = z * (LK + m * LB2) - m * LB2
    print("LKND =", LKND)
    
    PiKND_full = kurs.calc_PiKND(TBfull, LKND, effkLpcFull)["PiKND_full"]
    print("PiKND* =",PiKND_full)
    TKNDfull = kurs.calc_PiKND(TBfull, LKND, effkLpcFull)["TKNDfull"]
    print("TKND* =",TKNDfull) 
    
    PiKVD_full = PiK_full / PiKND_full
    print("PiKVD* =",PiKVD_full)
    
    LKVD_adiabatic = kurs.calc_LKVD(TKNDfull, TKfull, PiKVD_full)["LKVD_adiabatic"]
    print("LKVD adiabatic =",LKVD_adiabatic)
    
    effkHpcFull = LKVD_adiabatic / (LK - LKND)
    print("effk HPC* =", effkHpcFull)

#ПЕРВЫЙ КОНТУР

#Cечение КНД-КВД
    cpLPC = kurs.get_cp_real(TKNDfull)["cp"]
    kLPC = kurs.get_cp_real(TKNDfull)["k"]
    print(cpLPC, kLPC)
    PKNDfull = PBfull * PiKND_full
    print("PKND* =",PKNDfull)
    TKND = TKNDfull - cLPC *cLPC / 2 / cpLPC
    print("TKND =",TKND)
    PKND = PKNDfull * math.pow(TKNDfull / TKND, kLPC / (1 - kLPC))
    print("PKND =",PKND)
    roKND = calc_density(PKND, TKND, 287.0)
    print("roKND =",roKND)

#Сечение К-К
    cpHPC = kurs.get_cp_real(TKfull)["cp"]
    kHPC = kurs.get_cp_real(TKfull)["k"]
    print(cpHPC, kHPC)
    PKfull = PKNDfull * PiKVD_full
    print("PK* =",PKfull)
    LKVD = LKVD_adiabatic / effkHpcFull
    print("LKVD =",LKVD)
    TK = TKfull - cHPC *cHPC / 2 / cpHPC
    print("TK =",TK)
    PK = PKfull * math.pow(TKfull / TK, kHPC / (1 - kHPC))
    print("PK =",PK)
    roK = calc_density(PK, TK, 287.0)
    print("roK =",roK)
    
#Сечение Г-Г
    cpGas = kurs.get_cpmix(0.86, 0.14, TGfull, alpha)["cp_mix"]
    kGas = kurs.get_cpmix(0.86, 0.14, TGfull, alpha)["k"]
    RGas = kurs.get_cpmix(0.86, 0.14, TGfull, alpha)["R_mix"]
    print(cpGas, kGas)
    PGfull = PKfull * sigmaCC
    print("PG* =",PGfull)
    TG = TGfull - cG *cG / 2 / cpGas
    print("TG =",TG)
    PG = PGfull * math.pow(TGfull / TG, kGas / (1 - kGas))
    print("PG =",PG)
    roG = calc_density(PG, TG, 287.0)
    print("roG =",roG)
    #уточнение охлаждения, коррекция температуры Тг*
    TGcorr = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["TGcorr"]
    cpMix = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
    kMix = kurs.get_corr_parameters(TGfull, TKfull, cpGas, cpHPC, RGas, alpha, ksiTake, q_T, gAirBack)["k"]
    print(TGcorr, cpMix, kMix)
    
#Сечение ТВД-ТВД
    LTVD = LKVD / (1 + q_T + v_take)
    print("LTVD =",LTVD)
    TTVDfull = kurs.calc_LTurb(TGcorr, alpha, LTVD, effkHptFull)["TTurb"]
    print("TTVD* =",TTVDfull)
    PiTVD_full = kurs.calc_LTurb(TGcorr, alpha, LTVD, effkHptFull)["PiTurb"]
    print("PiTVD* =",PiTVD_full)
    PTVDfull = PGfull / PiTVD_full
    print("PTVD* =", PTVDfull)
    kMixTVDfull = kurs.get_cpmix(0.86, 0.14, TTVDfull, alpha)["k"]
    RMixTVDfull = kurs.get_cpmix(0.86, 0.14, TTVDfull, alpha)["R_mix"]
    print(kMixTVDfull, RMixTVDfull)
    TTVD = TTVDfull / (1 + (kMixTVDfull - 1) * MHPT * MHPT / 2)
    print("TTVD =", TTVD)
    PTVD = PTVDfull * math.pow(TTVDfull / TTVD, kMixTVDfull / (1 - kMixTVDfull))
    print("PTVD =", PTVD)
    kMixTVD = kurs.get_cpmix(0.86, 0.14, TTVD, alpha)["k"]
    RMixTVD = kurs.get_cpmix(0.86, 0.14, TTVD, alpha)["R_mix"]
    roTVD = calc_density(PTVD, TTVD, RMixTVD)
    print("roTVD =",roTVD)
    print(kMixTVD, RMixTVD)
    cTVD = MHPT * math.sqrt(kMixTVD * RMixTVD * TTVD)
    print("cTVD =", cTVD)
    
#Сечение ТНД-ТНД
    LTND = (LB2 * m + LKND) / (1 + q_T - v_take)
    print("LTND =", LTND)
    TTNDfull = kurs.calc_LTurb(TTVDfull, alpha, LTND, effkLptFull)["TTurb"]
    print("TTND* =",TTNDfull)
    PiTND_full = kurs.calc_LTurb(TTVDfull, alpha, LTND, effkLptFull)["PiTurb"]
    print("PiTND* =",PiTND_full)
    PTNDfull = PTVDfull / PiTND_full
    print("PTND* =", PTNDfull)
    kMixTNDfull = kurs.get_cpmix(0.86, 0.14, TTNDfull, alpha)["k"]
    RMixTNDfull = kurs.get_cpmix(0.86, 0.14, TTNDfull, alpha)["R_mix"]
    print(kMixTNDfull, RMixTNDfull)
    TTND = TTNDfull / (1 + (kMixTNDfull - 1) * MLPT * MLPT / 2)
    print("TTND =", TTND)
    PTND = PTNDfull * math.pow(TTNDfull / TTND, kMixTNDfull / (1 - kMixTNDfull))
    print("PTND =", PTND)
    kMixTND = kurs.get_cpmix(0.86, 0.14, TTND, alpha)["k"]
    RMixTND = kurs.get_cpmix(0.86, 0.14, TTND, alpha)["R_mix"]
    roTND = calc_density(PTND, TTND, RMixTND)
    print("roTND =",roTND)
    print(kMixTND, RMixTND)
    cTND = MLPT * math.sqrt(kMixTND * RMixTND * TTND)
    print("cTND =", cTND)
    
#Сечение С1-С1
    TC1full = TTNDfull
    print("TC1* =", TC1full)
    kMixC1 = kurs.get_cpmix(0.86, 0.14, TC1full, alpha)["k"]
    RMixC1 = kurs.get_cpmix(0.86, 0.14, TC1full, alpha)["R_mix"]
    cc1t = math.sqrt(2 * kMixC1 / (kMixC1 + 1) * RMixC1 * TC1full)
    print("cc1t =", cc1t)
    cc1 = phiC1 * cc1t
    print("cc1 =", cc1)
    sigmaC1 = kurs.GDF_pressure(kMixC1, 1, 1)["GDF_P"] / kurs.GDF_pressure(kMixC1, 1, phiC1)["GDF_P"]
    print("sigmaC1 =", sigmaC1)
    PiC1_full = PTNDfull * sigma1 / P0
    print("PiC1* =", PiC1_full)
    PiCr = math.pow((kMixC1 + 1) / 2, kMixC1 / (kMixC1 - 1))
    print("PiCr* =", PiCr)
    cc1_docr = kurs.calc_cc1_docr(TC1full, T0, PiC1_full, phiC1, alpha)["cc1_docr"]
    print("cc1 =", cc1_docr)
    PC1full = PTNDfull * sigmaC1
    print("PC1* =", PC1full)
    PC1 = P0
    print("PC1 =", PC1)
    cpMixC1ave = kurs.get_cpmix_ave(0.86, 0.14, TC1full, T0, alpha)["cp_mix"]
    RMixC1ave = kurs.get_cpmix_ave(0.86, 0.14, TC1full, T0, alpha)["R_mix"]
    TC1 = TC1full - cc1_docr * cc1_docr / 2 / cpMixC1ave
    print("TC1 =", TC1)
    roC1 = calc_density(PC1, TC1, RMixC1ave)
    print("roC1 =",roC1)
    
#ВТОРОЙ КОНТУР

#Сечение В2-В2
    PB2full = PBfull * 1.7
    print("PB2* =", PB2full)
    cpB2 = kurs.get_cp_real(TB2full)["cp"]
    kB2 = kurs.get_cp_real(TB2full)["k"]
    TB2 = TB2full - cB2 * cB2 / 2 / cpB2
    print("TB2 =", TB2)
    PB2 = PB2full * math.pow(TB2full / TB2, kB2 / (1 - kB2))
    print("PB2 =", PB2)
    roB2 = calc_density(PB2, TB2, 287.0)
    print("roB2 =", roB2)
    
#Сечение С2-С2
    TC2full = TB2full
    print("TC2* =", TC2full)
    PC2full = PB2full * sigma2
    print("PC2* =", PC2full)
    kC2 = kurs.get_cp_real(TC2full)["k"]
    sigmaC2 = kurs.GDF_pressure(kC2, 1, 1)["GDF_P"] / kurs.GDF_pressure(kC2, 1, phiC2)["GDF_P"]
    print("sigmaC2 =", sigmaC2) #дать комментарий почему так
    PiC2_full = PC2full * sigma2 / P0 / sigmaC2
    print("PiC2* =", PiC2_full)
    Pi2Cr = math.pow((kC2 + 1) / 2, kC2 / (kC2 - 1))
    print("Pi2Cr* =", Pi2Cr)
    cpC2 = kurs.get_cp_air(TC2full, T0)["cp"]
    kC2 = kurs.get_cp_air(TC2full, T0)["k"]
    cc2_docr = kurs.calc_cc2_docr(TC2full, T0, PiC2_full, phiC2)["cc2_docr"]
    print("cc2 =", cc2_docr)
    PC2 = P0
    print("PC2 =", PC2)
    TC2 = TC2full - cc2_docr * cc2_docr / 2 / cpC2
    print("TC2 =", TC2)
    roC2 = calc_density(PC2, TC2, 287.0)
    print("roC2 =", roC2)
    
#Основные параметры двигателя
    P_spec = cc1 * (1 + q_T + v_take) / (m + 1) + cc2_docr * m / (m + 1)
    print("P spec =", P_spec)
    Gair = 230000 * 1.07 / P_spec
    print("Gair =", Gair)
    Gair1 = Gair / (1 + m)
    print("Gair1 =", Gair1)
    Gair2 = Gair - Gair1
    print("Gair2 =", Gair2)
    Ggas = (1 + q_T + v_take) * Gair1
    print("Ggas =", Ggas)
    C_spec = (3600 * q_T * (1 - ksiTake)) / ((1 + m) * P_spec)
    print("C_spec =", C_spec)
    EngEff = ((1 + q_T - v_take) * (cc1_docr * cc1_docr / 2) + (m * cc2_docr * cc2_docr / 2)) / (Hu * q_T * effkGas)
    print("Engine eff =", EngEff)
    NB2 = LB2 * Gair2
    print("NB2 =", NB2)
    NKND = LKND * Gair1
    print("NKND =", NKND)
    NKVD = LKVD * Gair1
    print("NKVD =", NKVD)
    NTVD = LTVD * Ggas
    print("NTVD =", NTVD)
    NTND = LTND * Ggas
    print("NTND =", NTND)

    temp = [1 * 2,2 * 2,3 * 2,5 * 2,10 * 2,3 * 2,4 * 2,7 * 2,3 * 2]
    x = formArray(temp)
    labels = ['Н',"ВХ","В","КНД","К","Г","ТВД","ТНД","С1"]
    

# Построение второго графика на правой оси
    
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8*2.5, 6*2.5))
    first_axes = axes[0]
    first_axes.set_xticks(x)
    first_axes.set_xticklabels(labels, rotation=0, ha='center')
    y1 = [T0,TBXfull,TBfull,TKNDfull,TKfull,TGfull,TTVDfull,TTNDfull,TC1full]
    y2 = [T0,TBX,TB,TKND,TK,TG,TTVD,TTND,TC1]
    first_axes.plot(x, y1, 'red')
    first_axes.set_ylabel('Температуры')
    first_axes.tick_params(axis='y', labelcolor='b')
    first_axes.set_ylim(200,1800)
    right_ax = first_axes.twinx()
    right_ax.plot(x, y2, 'red',linestyle='-.')
    right_ax.set_ylim(200,1800)
    for x_i in x:
        plt.axvline(x=x_i, color='black', linestyle='--')
        
    second_axes = axes[1]
    y1 = [P0,PBXfull,PBfull,PKNDfull,PKfull,PGfull,PTVDfull,PTNDfull,PC1full]
    y2 = [P0,PBX,PB,PKND,PK,PG,PTVD,PTND,PC1]
    second_axes.set_xticks(x)
    second_axes.set_xticklabels(labels, rotation=0, ha='center')
    second_axes.plot(x, y1, 'blue')
    second_axes.set_ylabel('Давления')
    second_axes.tick_params(axis='y', labelcolor='b')
    second_axes.set_ylim(20*1000,2800*1000)
    right_ax_sec = second_axes.twinx()
    right_ax_sec.plot(x, y2, 'blue',linestyle='-.')
    right_ax_sec.set_ylim(20*1000,2800*1000)
    for x_i in x:
        plt.axvline(x=x_i, color='black', linestyle='--')
        
    thir_axes = axes[2]
    y2 = [0,cBX,cB,cLPC,cHPC,cG,cTVD,cTND,cc1]
    thir_axes.set_xticks(x)
    thir_axes.set_xticklabels(labels, rotation=0, ha='center')
    thir_axes.plot(x, y2, 'black')
    thir_axes.set_ylabel('Скорости')
    thir_axes.tick_params(axis='y', labelcolor='b')
    for x_i in x:
        thir_axes.axvline(x=x_i, color='black', linestyle='--')
    plt.tight_layout() # Автоматически корректирует расположение

    # 4. Отображение графика
    plt.savefig('sections.png', dpi=300, bbox_inches='tight')
def formArray(arr):
    res = []
    temp = 0
    for value in arr:
        res.append(temp + value)
        temp += value
    return res
    
# Вызов функции main
if __name__ == "__main__":
    main()