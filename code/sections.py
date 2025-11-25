import math
import lab1
import matplotlib.pyplot as plt
import kursach as kurs


#ПРЕДРАСЧЁТ

#КОЭФФИЦИЕНТЫ
coef = {
        "g_c": 0.86,
        "sigma_intake": 0.99,
        "sigma_cc": 0.955,
        "sigma_1": 0.99,
        "phi_c1": 0.97,
        "phi_c2": 0.97,
        "effk_gas": 0.995,
        "effk_comp_full": 0.81,
        "effk_fan_full": 0.87,
        "a": 0.03,
        "effk_hpt_full": 0.91,
        "effk_lpt_full": 0.92,
        "ksi_take": 0.1,
        "g_air_back": 0.07,
    }

gC           =  0.86
sigmaBX      = 0.99
sigmaCC      = 0.955
sigma1       = 0.99
sigma2 = 0.98
phiC1        = 0.97
phiC2        = 0.97
effkGas      = 0.995
effkCompFull = 0.81
effkLpcFull = 0.87
effkFanFull  = 0.87
aCoeff       = 0.03
effkHptFull  = 0.91
effkLptFull  = 0.92
ksiTake      = 0.1
gAirBack     = 0.07
z = 0.6

cB   = 220.0
cBX  = cB - 45.0 #45
cKND = 180.0
cB2  = cB - 17.5
cK = 155.0
cG   = 155.0 # 135
MHPT = 0.5 # 0.6
MLPT = 0.5 # 0.325
dB = (0.30 + 0.65) / 2
dHPC = (0.50 + 0.65) / 2

DcpDivHvyhHighPressure = (6.0 + 20.0) / 2
DcpDivHvyhOther = 3.0
F1           = (0.25 + 0.35) / 2
P0           = 101325.0
T0           = 288.0
ro0          = P0 / (287.0 * T0)

#ПРЕДРАСЧЁТ
m = 6.7 # 
TGfull = 1696.656
TKfull = 805.1776
PiK_full = 26.2500
variant = lab1.calc_proto(coef,TGfull, m,PiK_full)
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
cBX = 175.0
cB = 220.0
sigmaBX = 0.99

def main():
    
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
          "\n\troBX =", roBX,
        
    )
    #print(f"{TBX:.2f}|{PBX:.2f}|{PoBX:.2f}")

#Сечение В-В
    print("В-В:")
    PBfull = PBXfull * sigmaBX
    TBfull = TBXfull
    cP = 1004.5
    TB = TBfull - (cB * cB / cP)
    PB = PBfull * math.pow(TBfull/TB, (1.4/(1-1.4)))
    roB = calc_density(PB, TB, 287.0)
    print("\tTB* =", TBfull,
        "\n\tTB =", TB,
        "\n\tPB* =", PBfull,
        "\n\tPB =", PB,
        "\n\troB =", roB, 
    )
    
#Доп параметры для сечений
    print("Params:")
    compressor = kurs.calculate_compressor_temperature(
                288, PiK_full, coef["effk_comp_full"]
            )
    combustion_props = kurs.calculate_combustion_properties(
                0.86, 0.14, TKfull, TGfull, coef["effk_gas"]
            )
    q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
    print("\tq_T = ",q_T)
    Hu = combustion_props["Hu"]
    v_take = coef["ksi_take"] - coef["g_air_back"]
    LB2 = calculate_LB2(q_T,v_take,Xopt,Lfree,coef["effk_lpt_full"],6.7)
    print("\tLB2 = ",LB2)
    TB2full = kurs.calculate_fan_temperature(TBfull,LB2,coef["effk_fan_full"])[0]
    print("\tTB2full = ",TB2full)
    LB2 = kurs.calc_LB2(TBfull,1.7,coef["effk_fan_full"])["LB2"]
    TB2full = kurs.calc_LB2(TBfull,1.7,coef["effk_fan_full"])["TB2full"]
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
    print("\tkMixC1 = ",kMixC1)
    cc1t = math.sqrt(2 * kMixC1 / (kMixC1 + 1) * RMixC1 * TC1full)
    print("\tcc1t =", cc1t)
    cc1 = phiC1 * cc1t
    print("\tcc1 =", cc1)
    sigmaC1 = kurs.GDF_pressure(kMixC1, 1, 1)["GDF_P"] / kurs.GDF_pressure(kMixC1, 1, phiC1)["GDF_P"]
    GDF1 = kurs.GDF_pressure(kMixC1, 1, phiC1)["GDF_P"]
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
    kMixC1ave = kurs.get_cpmix_ave(TC1full, T0, alpha, ksiTake, q_T, gAirBack)["k"]
    cpMixC1ave = kurs.get_cpmix_ave(TC1full, T0, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
    RMixC1ave = kurs.get_cpmix_ave(TC1full, T0, alpha, ksiTake, q_T, gAirBack)["R_mix"]
    print("\tkMixC1ave = ",kMixC1ave,"\n\tcpMixC1ave = ", cpMixC1ave,"\n\tRMixC1ave = ", RMixC1ave)
    TC1 = TC1full - cc1_docr * cc1_docr / 2 / cpMixC1ave
    print("\tTC1 =", TC1)
    roC1 = calc_density(PC1, TC1, RMixC1ave)
    print("\troC1 =",roC1)
    
#ВТОРОЙ КОНТУР

#Сечение В2-В2
    print("В2-В2:")
    PB2full = PBfull * 1.7
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
    Gair = 230000 * 1.07 / P_spec
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
    second_axes.set_ylim(20*1000,2800*1000)
    right_ax_sec = second_axes.twinx()
    right_ax_sec.plot(x, y2, label = 'Статика',color = 'blue',linestyle='-.')
    right_ax_sec.set_ylim(20*1000,2800*1000)
    
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

    #Геометрия
    # Исходные параметры (заполни нужными значениями)
    F_BX = Gair / (cBX * roBX)
    # Диаметр круглого входного сечения:
    D_BX = math.sqrt(4 * F_BX / math.pi)
    
        # Исходные параметры
# Вентилятор
    var_calc_B = 0
    Dcp_B = 0.0
    Dvt_B = 0.0
    d_vt_B = 0.35
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
    #Внешний
    #
    # Площадь общего сечения перед разделителем (сечение P3J)
    F_rzd = Gair / (cB2 * roB2)
    # Постоянный средний диаметр
    D_cp_rzd = Dcp_B  # const
    # Диаметр корпуса и втулки в сечении разделителя
    D_vt_rzd = D_cp_rzd - (F_rzd / (math.pi * D_cp_rzd))
    D_V_rzd = 2 * D_cp_rzd - D_vt_rzd
    # Высота лопатки СА вентилятора (или мнимая высота)
    h_SA_rzd = (D_V_rzd - D_vt_rzd) / 2
    #внутренний
    #
    # Площадь сечения канала второго контура
    F_B2 = Gair2 / (cB2 * roB2)
    # Диаметр разделителя
    D_rzd = math.sqrt(D_V_rzd**2 - (4 * F_B2 / math.pi))
    # Высота канала второго контура или высота лопатки НА вентилятора
    h_B2 = (D_V_rzd - D_rzd) / 2
    # Толщина разделителя
    b_rzd = 0.085 * h_SA_rzd  # ... или 0.1 * h_SA_red
    # Внешний диаметр канала первого контура
    F1 = Gair1 / (cB2 * roB2)
    D_vknd = D_rzd - 2 * b_rzd
    # Втулочный диаметр канала первого контура
    D_vt_vknd = math.sqrt(D_vknd ** 2 - (4 * F1 / math.pi))
    # Средний диаметр канала первого контура
    D_cp_vknd = (D_vknd + D_vt_vknd) / 2
    # Высота лопатки на входе в канал первого контура
    h_vknd = (D_vknd - D_vt_vknd) / 2
    # Диаметр корпуса и втулки на выходе из КНД
#КНД
    # Площадь кольцевого сечения КНД
    F_knd = Gair1 / (cKND * roKND)    
    D_cp_knd = D_cp_vknd
    D_vt_knd = D_cp_knd - (F_knd / math.pi / D_cp_knd)
    D_knd = (2 * D_cp_knd - D_vt_knd)
    # Высота лопатки последнего СА КНД
    h_knd = (D_knd - D_vt_knd) / 2
    # Проверка допустимости геометрии КНД
    # 1. Проверка минимальной длины лопатки
    min_blade_length = 0.012  # 12 мм в метрах
    is_min_length_ok = h_knd >= min_blade_length
    # 2. Проверка относительного диаметра втулки
    relative_hub_diameter = D_vt_knd / D_knd
    max_relative_hub_diameter = 0.92  # максимально допустимое значение
    min_relative_hub_diameter = 0.82  # минимальное рекомендуемое значение
    is_hub_diameter_ok = (relative_hub_diameter <= 0.82)
    

    # Общая проверка геометрии КНД
    is_knd_geometry_valid = is_min_length_ok and is_hub_diameter_ok
##вКВД-вКВД
    F_vkvd = F_knd #Так как пи не равны
    d_vt_kvd = 0.6 #!!!!!!!!!!!!!НЕУВЕРЕН
    D_vkvd = math.sqrt((4 * F_vkvd) / (math.pi * (1 - d_vt_kvd * d_vt_kvd)))
    # Диаметр втулки КВД
    D_vt_vkvd = D_vkvd * d_vt_kvd
    # Средний диаметр КВД
    D_cp_vkvd = (D_vkvd + D_vt_vkvd) / 2
    h_vkvd = (D_vkvd - D_vt_vkvd) / 2
#K-K
    # Площадь сечения на выходе из КВД
    F_K = Gair1 / (cK * roK) #!!!!!!!!!! почему cK а не cK вот roK норм
    # Постоянный внешний диаметр канала
    D_K = D_vkvd  # const
    # Диаметр втулки на выходе из КВД
    D_vt_K = math.sqrt(D_K**2 - (4 * F_K / math.pi))
    # Средний диаметр на выходе из КВД
    D_cp_K = (D_K + D_vt_K) / 2
    h_K = (D_K - D_vt_K) / 2
    #Проверка
    compressor_blade_check = h_K >= 0.015
    rel_vt_vkvd_check = (D_vt_vkvd / D_vkvd) >= 0.5
    rel_vt_K_check = 0.92 >= (D_vt_K / D_K) >= 0.87
    rel_blade_check = 6.25 >= (h_vknd / h_K) >= 2
    compressor_check = (compressor_blade_check
                        and rel_vt_vkvd_check
                        and rel_vt_K_check
                        and rel_blade_check)
##ТВД-ТВД
    #Dср твд / h твд = 13
    D_cp_H = 11
    F_tvd = Gair1 * (1 + q_T - v_take) / (cTVD * roTVD)
    h_tvd = math.sqrt((4 * F_tvd) / (math.pi*(math.pow((D_cp_H + 1),2) - math.pow((D_cp_H - 1),2))))
    D_tvd = h_tvd * (D_cp_H + 1)
    Dvt_tvd = h_tvd * (D_cp_H - 1)
    Dcp_tvd = D_cp_H * h_tvd
    # Площадь проходного сечения Г-Г (без учета возврата воздуха)
    F_g = (Gair1 * (1 + q_T - ksiTake)) / (cG * roG)
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
    var_calc_tnd = 0 # 0-вт,1-ср,2-корпус
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
    Fc1 = (Gair1 * (1 + q_T - ksiTake)) / (cc1 * roC1)
    Dc1 = math.sqrt(4 * Fc1 / math.pi)
    #вС2 надо сначала посчитать все остальное
    #C2-C2
    Fc2 = Gair2 / (cc2_docr * roC2)
    ##Dc2 = math.sqrt((4 * Fc2) / math.pi + Dvt_c2)
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
        print(f"Средний диаметр D_cp_rzd: {D_cp_rzd:.4f} м")
        print(f"Диаметр корпуса D_V_rzd: {D_V_rzd:.4f} м")
        print(f"Диаметр втулки D_vt_rzd: {D_vt_rzd:.4f} м")
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
        print(f"Средний диаметр D_cp_knd: {D_cp_knd:.4f} м")
        print(f"Диаметр корпуса D_knd: {D_knd:.4f} м")
        print(f"Диаметр втулки D_vt_knd: {D_vt_knd:.4f} м")
        print(f"Высота лопатки h_knd: {h_knd:.4f} м")
        print(f"Относительный диаметр втулки: {relative_hub_diameter:.3f}")
        print(f"Геометрия КНД корректна: {'ДА' if is_knd_geometry_valid else 'НЕТ'}")
        
        print("\n" + "=" * 30)
        print("КВД - ВХОД (vKVD)")
        print("=" * 30)
        print(f"Относительный диаметр втулки d_vt_kvd: {d_vt_kvd:.3f}")
        print(f"Внешний диаметр D_vkvd: {D_vkvd:.4f} м")
        print(f"Диаметр втулки D_vt_vkvd: {D_vt_vkvd:.4f} м")
        print(f"Средний диаметр D_cp_vkvd: {D_cp_vkvd:.4f} м")
        print(f"Высота лопатки h_vkvd: {h_vkvd:.4f} м")

        print("\n" + "=" * 30)
        print("КВД - ВЫХОД (K)")
        print("=" * 30)
        print(f"Площадь сечения F_K: {F_K:.4f} м²")
        print(f"Внешний диаметр D_K: {D_K:.4f} м")
        print(f"Диаметр втулки D_vt_K: {D_vt_K:.4f} м")
        print(f"Средний диаметр D_cp_K: {D_cp_K:.4f} м")
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
        print(f"Площадь сечения Fc1: {Fc1:.4f} м²")
        print(f"Диаметр сопла Dc1: {Dc1:.4f} м")

        print("\n" + "=" * 50)
        print("ПРОВЕРКИ ГЕОМЕТРИИ")
        print("=" * 50)
        print(f"КНД - мин. высота лопатки ({min_blade_length*1000:.0f} мм): {'СОБЛЮДЕНО' if is_min_length_ok else 'НЕ СОБЛЮДЕНО'}")
        print(f"КНД - отн. диаметр втулки ({min_relative_hub_diameter}-{max_relative_hub_diameter}): {'СОБЛЮДЕНО' if is_hub_diameter_ok else 'НЕ СОБЛЮДЕНО' }")
        if compressor_check:
            print(f"ГЕОМЕТРИЯ КОМПРЕССОРА: СОБЛЮДЕНА")
        else:
            print(f"КОМПРЕССОР ЧЕК:\nВысота лопатки последнего СА:{"СОБЛЮДЕНО" if compressor_blade_check else "НЕ СОБЛЮДЕНО"}")
            print(f"Относительный диаметр втулки в сечении К:{"СОБЛЮДЕНО" if rel_vt_K_check else "НЕ СОБЛЮДЕНО - "}{D_vt_vkvd / D_vkvd}")
            print(f"Относительный диаметр втулки в сечении вКВД:{"СОБЛЮДЕНО" if rel_vt_vkvd_check else "НЕ СОБЛЮДЕНО"}")
            print(f"отношение длин лопаток:{"СОБЛЮДЕНО" if rel_blade_check else "НЕ СОБЛЮДЕНО"}")
        if turbine_check:
            print(f"ГЕОМЕТРИЯ ТУРБИНЫ: СОБЛЮДЕНА")
        else:
            print(f"ТУРБИНА ЧЕК:\nОтношение длин лопаток ТВД:{"СОБЛЮДЕНО" if tvd_blade_check else "НЕ СОБЛЮДЕНО"}")
            print(f"Отношение длин лопаток ТНД:{"СОБЛЮДЕНО" if rel_blade_tnd else "НЕ СОБЛЮДЕНО"}")
            print(f"Приведенная длина лопатки в сечении ТНД:{"СОБЛЮДЕНО" if rel_blade_check else "НЕ СОБЛЮДЕНО"}")
            
    printGeometry()
    figNew, axesNew = plt.subplots(nrows=1, ncols=1, figsize=(8*2.5, 6*2.5))
    first_axes_new = axesNew
    temp = [1,2,3,6,7,10,12,13,15]
    x = temp
    labels = ['В',"Рзд","КНДв","КНДвых","КВДвх","КВДвых","Г","ТВД","ТНД"]
    first_axes_new.set_xticks(x)
    first_axes_new.set_xticklabels(labels, rotation=0, ha='center')
    Dvt = [Dvt_B,D_vt_rzd,D_vt_vknd,D_vt_knd,D_vt_vkvd,D_vt_K,Dvt_g,Dvt_tvd,Dvt_tnd]
    Dcp = [Dcp_B,D_cp_rzd,D_cp_vknd,D_cp_knd,D_cp_vkvd,D_cp_K,Dcp_g,Dcp_tvd,Dcp_tnd]
    D = [D_B,D_V_rzd,D_vknd,D_knd,D_vkvd,D_K,D_g,D_tvd,D_tnd]
    first_axes_new.plot(x, Dvt, label='Втулка', color='black')
    first_axes_new.plot(x, Dcp, label='Средняя', color='orange', linestyle='-.')
    first_axes_new.plot(x, D, label='Корпус', color='black')
    first_axes_new.legend(loc='best')
    for x_i in x:
        plt.axvline(x=x_i,color='black', linestyle='--')
    plt.savefig("../diagrams/comp.png")
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