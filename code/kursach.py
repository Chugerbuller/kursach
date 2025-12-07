import numpy as np
import math


def get_cp_air(T1, T2):
    poly =[-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]

    # Вычисление полинома (коэффициенты идут от старшей степени к младшей)
    return get_cp(T1,T2,poly)

def get_cp(T1, T2, poly):
    """
    Вычисляет среднюю теплоемкость воздуха в заданном температурном диапазоне.
    
    Параметры:
    T1 - начальная температура [K]
    T2 - конечная температура [K]
    
    Возвращает:
    cp - средняя теплоемкость воздуха [Дж/(кг·К)]
    """
    R = 287.0
    n = len(poly)
    cp = 0
    
    for i in range(n):
        power = n - i - 1
        cp += poly[i] * (T2**(power + 1) - T1**(power + 1)) / (power + 1)
    
    cp = cp / (T2 - T1)
    cv = cp - R
    k = cp/cv
    return {"cp": cp, "cv": cv, "k": k}

def get_cp_component(T1, T2, poly, R):
    """
    Вычисляет среднюю теплоемкость воздуха в заданном температурном диапазоне.
    
    Параметры:
    T1 - начальная температура [K]
    T2 - конечная температура [K]
    
    Возвращает:
    cp - средняя теплоемкость воздуха [Дж/(кг·К)]
    """
    n = len(poly)
    cp = 0
    
    for i in range(n):
        power = n - i - 1
        cp += poly[i] * (T2**(power + 1) - T1**(power + 1)) / (power + 1)
    
    cp = cp / (T2 - T1)
    cv = cp - R
    k = cp/cv
    return {"cp": cp, "cv": cv, "k": k}


def get_cp_real(T):
    poly =[-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]

    # Вычисление полинома (коэффициенты идут от старшей степени к младшей)
    return get_cp2(T,poly)

def get_cp2(T, poly):
    """
    Вычисляет истинную теплоемкость воздуха заданной температуры.
    
    Параметры:
    T - температура [K]
    
    Возвращает:
    cp - истинная теплоемкость воздуха [Дж/(кг·К)]
    """
    R=287.0
    n = len(poly)
    cp = 0
    
    for i in range(n):
        power = n - i - 1
        cp += poly[i] * T ** power
    
    cp = cp
    cv = cp - R
    k = cp/cv
    return {"cp": cp, "cv": cv, "k": k}

def get_cp3(T, poly, R):
    """
    Вычисляет истинную теплоемкость компонента заданной температуры.
    
    Параметры:
    T - температура [K]
    
    Возвращает:
    cp - истинная теплоемкость компонента [Дж/(кг·К)]
    """

    n = len(poly)
    cp = 0
    
    for i in range(n):
        power = n - i - 1
        cp += poly[i] * T ** power
    
    cp = cp
    cv = cp - R
    k = cp/cv
    return {"cp": cp, "cv": cv, "k": k}


def calculate_compressor_temperature(TH, pik, effK, R=287.0, max_iter=1000, tol=1e-16): #норм
    if TH <= 0 or pik <= 0 or effK <= 0:
        raise ValueError("Параметры должны быть положительными")

    TK_old = 0.0
    TK = 300.0
    iter_count = 0
    converged = False

    # История итераций для отладки
    history = []

    while abs(TK - TK_old) > tol and iter_count < max_iter:
        try:
            cpAir = get_cp_air(TH,TK)

            TK_old = TK
            TK = TH * (1 + (pik ** ((cpAir['k'] - 1) / cpAir['k']) - 1) / effK)

            # Сохраняем историю итераций
            history.append(
                {"iteration": iter_count, "TK": TK, "cpAir": cpAir['cp'], "kAir": cpAir['k']}
            )

            iter_count += 1

        except (ValueError, ZeroDivisionError) as e:
            print(f"Ошибка на итерации {iter_count}: {e}")
            break

    converged = abs(TK - TK_old) <= tol

    return {
        "TK": TK,
        "cpAir": cpAir['cp'],
        "cvAir": cpAir['cv'],
        "kAir": cpAir['k'],
        "iterations": iter_count,
        "converged": converged,
        "history": history,
    }


def calculate_combustion_properties(gC, gH, TK, TG, effG, max_iter=1000, tol=1e-16):#норм
    """
    Вычисление свойств продуктов сгорания с итерационным уточнением

    Parameters:
    C : float
        Доля углерода в топливе [кг/кг]
    H : float
        Доля водорода в топливе [кг/кг]
    TK : float
        Температура после компрессора [K]
    TG : float
        Температура горения [K]
    effG : float
        Эффективность сгорания топлива
    max_iter : int, optional
        Максимальное количество итераций
    tol : float, optional
        Допустимая погрешность сходимости

    Returns:
    dict : результаты расчета
    """
    # 1. Вычисляем удельные газовые постоянные
    R_CO2 = 8314.2 / (12 + 16 * 2)  # для углекислого газа
    R_H2O = 8314.2 / (2 + 16)  # для водяного пара
    R_N2 = 8314.2 / 28  # для азота
    R_O2 = 8314.2 / 32  # для кислорода

    # 2. Вычисляем низшую теплоту сгорания и теоретическое количество воздуха
    Hu = (33800 * gC + 102500 * gH) * 1000 # Дж/кг
    L0 = (8 / 3 * gC + 8 * gH) / 0.23  # кг/кг

    # 3. Функции для расчета теплоемкостей компонентов
    # (аналоги get_cp_air для каждого компонента)

    # 4. Итерационный процесс
    k_old = 1.4  # начальное предположение
    k = 0.0
    alpha = 1.0  # начальный коэффициент избытка воздуха
    iter_count = 0
    converged = False
    cp_CO2 = get_cp_component(TK, TG, [-5.2735e-11, 3.9194e-7, -1.1213e-3, 1.5466, 471.75], R_CO2)['cp']
    cp_H2O = get_cp_component(TK, TG, [8.2542e-11, -5.3927e-7, 1.0936e-3, -1.9361e-1,1842.53], R_H2O)['cp']
    cp_N2 = get_cp_component(TK, TG, [-3.5780e-14, 2.9022e-10, -8.8233e-7, 1.1757e-3,-4.7731e-1,1095.68], R_N2)['cp']
    cp_O2 = get_cp_component(TK, TG, [-4.7303e-14, 3.3563e-10, -8.4931e-7, 8.5606e-4, -1.0201e-1, 897], R_O2)['cp']
    
    #cp_CO2 = get_cp(TK, TG, [-5.2735e-11, 3.9194e-7, -1.1213e-3, 1.5466, 471.75])['cp']
    #cp_H2O = get_cp(TK, TG, [8.2542e-11, -5.3927e-7, 1.0936e-3, -1.9361e-1,1842.53])['cp']
    #cp_N2 = get_cp(TK, TG, [-3.5780e-14, 2.9022e-10, -8.8233e-7, 1.1757e-3,-4.7731e-1,1095.68])['cp']
    #cp_O2 = get_cp(TK, TG, [-4.7303e-14, 3.3563e-10, -8.4931e-7, 8.5606e-4, -1.0201e-1, 897])['cp']
    
    while abs(k - k_old) > tol and iter_count < max_iter:
        # Вычисляем массовые доли компонентов
        g_CO2 = (11 * gC) /(3 *(1 + alpha * L0))
        g_H2O = (9 * gH) / (1 + alpha * L0)
        g_N2 = (0.77 * L0 * alpha) / (1 + L0 * alpha)
        g_O2 = (0.23 * (alpha - 1) * L0) / (1 + L0 * alpha)
        # Вычисляем средние теплоемкости в интервале [TK, TG]
        history = []
        # cp_air = get_cp(TK, TG, np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]))

        # Вычисляем теплоемкость и газовую постоянную смеси
        cp_mix = cp_CO2 * g_CO2 + cp_H2O * g_H2O + cp_N2 * g_N2 + cp_O2 * g_O2

        R_mix = R_CO2 * g_CO2 + R_H2O * g_H2O + R_N2 * g_N2 + R_O2 * g_O2

        cv_mix = cp_mix - R_mix

        # Обновляем показатель адиабаты
        k_old = k
        k = cp_mix / cv_mix

        # Обновляем коэффициент избытка воздуха
        alpha = 1 / L0 * (Hu * effG / cp_mix / (TG - TK) - 1)

        # Сохраняем историю
        history.append(
            {
                "iteration": iter_count,
                "alpha": alpha,
                "k": k,
                "cp_mix": cp_mix,
                "R_mix": R_mix,
                "g_CO2": g_CO2,
                "g_H2O": g_H2O,
                "g_N2": g_N2,
                "g_O2": g_O2,
            }
        )

        iter_count += 1

    # Проверяем сходимость
    converged = abs(k - k_old) <= tol

    return {
        "alpha": alpha,
        "k": k, 
        "cp_mix": cp_mix,
        "R_mix": R_mix,
        "cv_mix": cv_mix,
        "mass_fractions": {"CO2": g_CO2, "H2O": g_H2O, "N2": g_N2, "O2": g_O2},
        "iterations": iter_count,
        "converged": converged,
        "history": history,
        "Hu": Hu,
        "L0": L0,
    }
    
def get_cpgas(gC, gH, TG, alpha, max_iter=1000, tol=1e-16):#норм
    
    # 1. Вычисляем удельные газовые постоянные
    R_CO2 = 8314.2 / (12 + 16 * 2)  # для углекислого газа
    R_H2O = 8314.2 / (2 + 16)  # для водяного пара
    R_N2 = 8314.2 / 28  # для азота
    R_O2 = 8314.2 / 32  # для кислорода

    L0 = (8 / 3 * gC + 8 * gH) / 0.23  # кг/кг

    # 2. Итерационный процесс
    k_old = 0.0  # начальное предположение
    k = 1.4
    iter_count = 0
    converged = False
    cp_CO2 = get_cp3(TG, [-5.2735e-11, 3.9194e-7, -1.1213e-3, 1.5466, 471.75], R_CO2)['cp']
    cp_H2O = get_cp3(TG, [8.2542e-11, -5.3927e-7, 1.0936e-3, -1.9361e-1,1842.53], R_H2O)['cp']
    cp_N2 = get_cp3(TG, [-3.5780e-14, 2.9022e-10, -8.8233e-7, 1.1757e-3,-4.7731e-1,1095.68], R_N2)['cp']
    cp_O2 = get_cp3(TG, [-4.7303e-14, 3.3563e-10, -8.4931e-7, 8.5606e-4, -1.0201e-1, 897], R_O2)['cp']
    
    while abs(k - k_old) > tol and iter_count < max_iter:
        # Вычисляем массовые доли компонентов
        g_CO2 = (11 * gC) /(3 *(1 + alpha * L0))
        g_H2O = (9 * gH) / (1 + alpha * L0)
        g_N2 = (0.77 * L0 * alpha) / (1 + L0 * alpha)
        g_O2 = (0.23 * (alpha - 1) * L0) / (1 + L0 * alpha)
        # Вычисляем средние теплоемкости в интервале [TK, TG]
        history = []
        # cp_air = get_cp(TK, TG, np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]))

        # Вычисляем теплоемкость и газовую постоянную смеси
        cp_mix = cp_CO2 * g_CO2 + cp_H2O * g_H2O + cp_N2 * g_N2 + cp_O2 * g_O2

        R_mix = R_CO2 * g_CO2 + R_H2O * g_H2O + R_N2 * g_N2 + R_O2 * g_O2

        cv_mix = cp_mix - R_mix

        # Обновляем показатель адиабаты
        k_old = k
        k = cp_mix / cv_mix

        # Сохраняем историю
        history.append(
            {
                "iteration": iter_count,
                "k": k,
                "cp_mix": cp_mix,
                "R_mix": R_mix,
                "g_CO2": g_CO2,
                "g_H2O": g_H2O,
                "g_N2": g_N2,
                "g_O2": g_O2,
            }
        )

        iter_count += 1

    # Проверяем сходимость
    converged = abs(k - k_old) <= tol

    return {
        "k": k, 
        "cp_mix": cp_mix,
        "R_mix": R_mix,
        "cv_mix": cv_mix,
        "mass_fractions": {"CO2": g_CO2, "H2O": g_H2O, "N2": g_N2, "O2": g_O2},
        "iterations": iter_count,
        "converged": converged,
        "history": history,
    }
    
def get_cpgas_ave(gC, gH, T1, T2, alpha, max_iter=1000, tol=1e-16):#норм
    
    # 1. Вычисляем удельные газовые постоянные
    R_CO2 = 8314.2 / (12 + 16 * 2)  # для углекислого газа
    R_H2O = 8314.2 / (2 + 16)  # для водяного пара
    R_N2 = 8314.2 / 28  # для азота
    R_O2 = 8314.2 / 32  # для кислорода

    L0 = (8 / 3 * gC + 8 * gH) / 0.23  # кг/кг

    # 2. Итерационный процесс
    k_old = 1.4  # начальное предположение
    k = 0.0
    iter_count = 0
    converged = False
    cp_CO2 = get_cp(T1, T2, [-5.2735e-11, 3.9194e-7, -1.1213e-3, 1.5466, 471.75])['cp']
    cp_H2O = get_cp(T1, T2, [8.2542e-11, -5.3927e-7, 1.0936e-3, -1.9361e-1,1842.53])['cp']
    cp_N2 = get_cp(T1, T2, [-3.5780e-14, 2.9022e-10, -8.8233e-7, 1.1757e-3,-4.7731e-1,1095.68])['cp']
    cp_O2 = get_cp(T1, T2, [-4.7303e-14, 3.3563e-10, -8.4931e-7, 8.5606e-4, -1.0201e-1, 897])['cp']
    
    while abs(k - k_old) > tol and iter_count < max_iter:
        # Вычисляем массовые доли компонентов
        g_CO2 = (11 * gC) /(3 *(1 + alpha * L0))
        g_H2O = (9 * gH) / (1 + alpha * L0)
        g_N2 = (0.77 * L0 * alpha) / (1 + L0 * alpha)
        g_O2 = (0.23 * (alpha - 1) * L0) / (1 + L0 * alpha)
        # Вычисляем средние теплоемкости в интервале [TK, TG]
        history = []
        # cp_air = get_cp(TK, TG, np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]))

        # Вычисляем теплоемкость и газовую постоянную смеси
        cp_mix = cp_CO2 * g_CO2 + cp_H2O * g_H2O + cp_N2 * g_N2 + cp_O2 * g_O2

        R_mix = R_CO2 * g_CO2 + R_H2O * g_H2O + R_N2 * g_N2 + R_O2 * g_O2

        cv_mix = cp_mix - R_mix

        # Обновляем показатель адиабаты
        k_old = k
        k = cp_mix / cv_mix

        # Сохраняем историю
        history.append(
            {
                "iteration": iter_count,
                "k": k,
                "cp_mix": cp_mix,
                "R_mix": R_mix,
                "g_CO2": g_CO2,
                "g_H2O": g_H2O,
                "g_N2": g_N2,
                "g_O2": g_O2,
            }
        )

        iter_count += 1

    # Проверяем сходимость
    converged = abs(k - k_old) <= tol

    return {
        "k": k, 
        "cp_mix": cp_mix,
        "R_mix": R_mix,
        "cv_mix": cv_mix,
        "mass_fractions": {"CO2": g_CO2, "H2O": g_H2O, "N2": g_N2, "O2": g_O2},
        "iterations": iter_count,
        "converged": converged,
        "history": history,
    }
    
def get_cpmix(T, alpha, ksiTake, q_T, gAirBack):
    gC = 0.86
    gH = 0.14
    g1 = (1 - ksiTake + q_T)
    g2 = gAirBack
    cpGas = get_cpgas(gC, gH, T, alpha, max_iter=1000, tol=1e-16)["cp_mix"]
    cpAir = get_cp_real(T)["cp"]
    RGas = get_cpgas(gC, gH, T, alpha, max_iter=1000, tol=1e-16)["R_mix"]
    cpMix = (cpGas * g1 + cpAir * g2) / (g1 + g2)
    RMix = (RGas * g1 + 287 * g2) / (g1 + g2)
    cvMix = cpMix - RMix
    kMix = cpMix/cvMix
    
    return {
        "k": kMix, 
        "cp_mix": cpMix,
        "R_mix": RMix,
        "cv_mix": cvMix,
    }
    
def get_cpmix_ave(T1, T2, alpha, ksiTake, q_T, gAirBack):
    gC = 0.86
    gH = 0.14
    g1 = (1 - ksiTake + q_T)
    g2 = gAirBack
    cpGas = get_cpgas_ave(gC, gH, T1, T2, alpha)["cp_mix"]
    cpAir = get_cp_air(T1, T2)["cp"]
    RGas = get_cpgas(gC, gH, T1, T2, alpha)["R_mix"]
    cpMix = (cpGas * g1 + cpAir * g2) / (g1 + g2)
    RMix = (RGas * g1 + 287 * g2) / (g1 + g2)
    cvMix = cpMix - RMix
    kMix = cpMix/cvMix
    
    return {
        "k": kMix, 
        "cp_mix": cpMix,
        "R_mix": RMix,
        "cv_mix": cvMix,
    }


def get_corr_parameters(TG, TK, cpGasG, cpAirK, RGas, alpha, ksiTake, q_T, gAirBack):
    TGcorr = TG
    TGcorrold = 0.0
    iter_count = 0
    gC = 0.86
    gH = 0.14
    g1 = (1 - ksiTake + q_T)
    g2 = gAirBack
    
    while abs(TGcorr - TGcorrold) > 1e-16 and iter_count < 1000:
        cpGas = get_cpgas(gC, gH, TGcorr, alpha, max_iter=1000, tol=1e-16)["cp_mix"]
        cpAir = get_cp_real(TGcorr)["cp"]
        cpMix = (cpGas * g1 + cpAir * g2) / (g1 + g2)
        RMix = (RGas * g1 + 287 * g2) / (g1 + g2)
        cvMix = cpMix - RMix
        kMix = cpMix/cvMix
        TGcorr = (cpGasG * g1 * TG + cpAirK * g2 * TK) / (cpGas * g1 + cpAir * g2)
        iter_count += 1
        
    return {
        "k": kMix, 
        "cp_mix": cpMix,
        "R_mix": RMix,
        "cv_mix": cvMix,
        "TGcorr": TGcorr
    }

def calculate_fan_temperature(TB, LB, effB, R=287.0):
    """Расчет температуры за вентилятором"""
    TB2old = 0.0
    TB2 = 300.0
    iter_count = 0

    while abs(TB2 - TB2old) > 1e-16 and iter_count < 100:
        cpAir = get_cp_air(TB, TB2)['cp']
        cvAir = cpAir - R
        kAir = cpAir / cvAir

        TB2old = TB2
        piB =  math.pow((effB * LB / cpAir / TB + 1),(kAir / (kAir - 1)))
        if piB > 1.95:
            piB = 1.95
        TB2 = TB * (1 + (math.pow(piB,((kAir - 1) / kAir)) - 1) / effB)

        iter_count += 1

    return TB2, piB, kAir, cpAir, iter_count

def calc_LB2(TB, piB_star, etaB_star):
    TB2old = 0.0
    TB2 = 300.0
    iter_count = 0
    
    while abs(TB2 - TB2old) > 1e-16 and iter_count < 100:
        cpAir = get_cp_air(TB, TB2)['cp']
        cvAir = cpAir - 287
        kAir = cpAir / cvAir

        TB2old = TB2
        LB2 =  cpAir * TB * (math.pow(piB_star,(kAir-1)/kAir) - 1) / etaB_star
        TB2 = TB * (1 + (math.pow(piB_star, (kAir-1)/kAir) - 1) / etaB_star)
        iter_count += 1
        
    return {"LB2":LB2, "TB2full" :TB2}

def calc_LK(TB, TK, piK, etaK):
    cpAir = get_cp_air(TB, TK)['cp']
    cvAir = cpAir - 287
    kAir = cpAir / cvAir
    
    LK = cpAir * TB * (math.pow(piK,(kAir-1)/kAir) - 1) / etaK
    
    return {"LK":LK}

def calc_LKND(z, LK, LB2, m, etaK, etaB, etaKND):
    LKND = (z * (LK * etaK + m * LB2 * etaB) - m * LB2 * etaB) / etaKND
    
    return {"LKND":LKND}

def calc_PiKND(TB, LKND, etaKND):
    TKNDold = 0.0
    TKND = 300.0
    iter_count = 0
    
    while abs(TKND - TKNDold) > 1e-16 and iter_count < 100:
        cpAir = get_cp_air(TB, TKND)['cp']
        cvAir = cpAir - 287
        kAir = cpAir / cvAir
        
        TKNDold = TKND
        PiKND_base = (LKND * etaKND /(cpAir * TB)) + 1
        PiKND = math.pow(PiKND_base, (kAir/(kAir-1)))
        TKND = TB * (1 + (math.pow(PiKND, (kAir-1)/kAir) - 1) / etaKND)
        iter_count += 1
        
    return {"PiKND_full":PiKND, "TKNDfull" :TKND}

def calc_LKVD(TKND, TK, piKVD):
    cpAir = get_cp_air(TKND, TK)['cp']
    cvAir = cpAir - 287
    kAir = cpAir / cvAir
    
    LKVD_adiabatic = cpAir * TKND * (math.pow(piKVD,(kAir-1)/kAir) - 1)
    return {"LKVD_adiabatic":LKVD_adiabatic}

def calc_LTurb(T1, alpha, LTurb, etaTurb, ksiTake, q_T, gAirBack):
    TTurbold = 0.0
    TTurb = 1500.0
    PiTurb = 5
    iter_count = 0
    
    while abs(TTurb - TTurbold) > 1e-16 and iter_count < 100:
        TTurbold = TTurb
        cp = get_cpmix_ave(T1, TTurb, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
        k = get_cpmix_ave(T1, TTurb, alpha, ksiTake, q_T, gAirBack)["k"]
        parent = (1 - math.pow(PiTurb, (1 - k) / k) )
        TTurb = T1 * (1 - parent * etaTurb)
        PiTurb = math.pow((1 - LTurb / (etaTurb * cp * T1)), k / (1 - k))
        iter_count += 1

    return{
        "TTurb": TTurb,
        "PiTurb": PiTurb
    }
    
def GDF_pressure(k, l, phi):
    parent = 1 - phi * phi * l * l * (k - 1)/(k + 1)
    GDF_P = math.pow(parent, k / (k - 1))
    return {"GDF_P": GDF_P}

def calc_cc1_docr(TC1, T0, Pi, Phi, alpha, ksiTake, q_T, gAirBack):
    cp = get_cpmix_ave(T0, TC1, alpha, ksiTake, q_T, gAirBack)["cp_mix"]
    k = get_cpmix_ave(T0, TC1, alpha, ksiTake, q_T, gAirBack)["k"]
    cc1_docr = Phi * math.sqrt(2 * cp * TC1 * (1 - math.pow(Pi, (1 - k) / k)))
    
    return{"cc1_docr": cc1_docr}

def calc_cc2_docr(TC2, T0, Pi, Phi):
    cp = get_cp_air(T0, TC2)["cp"]
    k = get_cp_air(T0, TC2)["k"]
    cc2_docr = Phi * math.sqrt(2 * cp * TC2 * (1 - math.pow(Pi, (1 - k) / k)))
    
    return{"cc2_docr": cc2_docr}
    
def get_in_geometry(DIN, F, geom_type):
    """
    Вычисляет геометрические параметры на основе относительного диаметра и площади сечения.
    """

    def case_0():
        D = math.sqrt(4 * F / (math.pi * (1 - DIN**2)))
        d = D * DIN
        DCP = (D + d) / 2
        return D, d, DCP

    def case_1():
        d = math.sqrt(4 * F / (math.pi * (1 / DIN**2 - 1)))
        D = d / DIN
        DCP = (D + d) / 2
        return D, d, DCP

    def case_2():
        DCP = math.sqrt(F * (1 + DIN) / (math.pi * (1 - DIN)))
        D = 2 * DCP / (1 + DIN)
        d = D * DIN
        return D, d, DCP

    switch = {0: case_0, 1: case_1, 2: case_2}

    if geom_type in switch:
        return switch[geom_type]()
    else:
        raise ValueError("Неверный тип геометрии. Допустимые значения: 0, 1, 2")
    import math


def get_out_geometry(DIN, F, geom_type):
    """
    Вычисляет геометрические параметры на основе определяющего диаметра и площади сечения.

    Параметры:
    DIN – определяющий диаметр (зависит от выбранной геометрии узла во входном сечении)
    F – площадь сечения
    geom_type – тип геометрии:
        0 - постоянный внешний диаметр
        1 - постоянный внутренний диаметр
        2 - постоянный средний диаметр

    Возвращает:
    D – внешний диаметр
    d – внутренний диаметр
    DCP – средний диаметр
    """

    match geom_type:
        case 0:
            # постоянный внешний диаметр. DIN – диаметр корпуса по входному сечению
            D = DIN
            d = math.sqrt(D**2 - 4 * F / math.pi)
            DCP = (D + d) / 2

        case 1:
            # постоянный внутренний диаметр. DIN – диаметр втулки по входному сечению
            d = DIN
            D = math.sqrt(d**2 + 4 * F / math.pi)
            DCP = (D + d) / 2

        case 2:
            # постоянный средний диаметр. DIN – средний диаметр по входному сечению
            DCP = DIN
            d = DCP - F / (math.pi * DCP)
            D = 2 * DCP - d

        case _:
            raise ValueError("Неверный тип геометрии. Допустимые значения: 0, 1, 2")

    return D, d, DCP


def calculate_turbine_stages(LTBD, effTBD, DCPTBD, DCPG, uK1, DBKBD):
    """
    Функция для расчета оптимального числа ступеней турбины на основе параметра Парсонса.

    Параметры:
    LTBD - работа турбины
    effTBD - КПД турбины
    DCPTBD - средний диаметр на входе в турбину
    DCPG - средний диаметр на выходе из турбины
    uK1 - скорость на входе в КВД на внешнем диаметре РК
    DBKBD - внешний диаметр на входе в КВД

    Возвращает:
    z - оптимальное число ступеней
    y - итоговое значение параметра Парсонса
    iter - количество выполненных итераций
    """

    # Теоретическая скорость истечения из турбины
    c0 = math.sqrt(LTBD / effTBD * 2)

    z = 1  # начальное число ступеней
    deltaD = (
        DCPTBD - DCPG
    )  # различие в средних диаметрах между входным и выходным сечениями турбины
    stepD = 0  # инициализация переменной шага по диаметру
    y = (
        uK1 / c0
    )  # начальный параметр Парсонса (скорость на среднем диаметре / теоретическая скорость)

    iter = 0  # счетчик итераций
    opt = True  # логический критерий остановки

    # Цикл подбора ступеней, пока не выполнится условие по параметру Парсонса
    while opt and iter < 100:
        # Проверяем, лежит ли критерий в заданном интервале, и меняем число ступеней
        if y >= 0.75 and (z - 1) > 0:
            z = z - 1
        if y <= 0.45:
            z = z + 1

        y = 0
        stepD = deltaD / z  # пересчитываем шаг между диаметрами

        # Вычисляем новые скорости по ступеням ТВД через отношение диаметров
        for it in range(1, z + 1):
            # Вычисляем скорость на текущей ступени через отношение диаметров
            uCPz = uK1 * (DCPG + stepD * it) / DBKBD
            # Суммируем отношения квадратов скоростей на среднем диаметре
            y += (uCPz / c0) ** 2

        y = math.sqrt(y)  # извлекаем корень и получаем итоговое значение Парсонса

        # Если итоговое значение легло в допустимый диапазон, число ступеней успешно подобрано
        if 0.45 < y < 0.75:
            opt = False

        iter += 1

    return z, y, iter



def optimal_parameter(phi_c1, phi_c2, q_T, v_take, m, eta_LPT_star, eta_F_star):#норм
    """
    Вычисляет оптимальный параметр x_opt по формуле:

    x_opt = 1 / [1 + (φ_c1² * (1 + q_T - v_take)) / (φ_c2² * m * η_THA* * η_B*)]

    Parameters:
    phi_c1 (float): Коэффициент скорости на входе (φ_c1)
    phi_c2 (float): Коэффициент скорости на выходе (φ_c2)
    q_T (float): Параметр q_T
    v_take (float): Параметр v_take
    m (float): Параметр m
    eta_LPT_star (float): КПД турбины низкого давления (η_THA*)
    eta_F_star (float): КПД вентилятора (η_B*)

    Returns:
    float: Оптимальный параметр x_opt
    """
    numerator = math.pow(phi_c1,2) * (1 + q_T - v_take)
    denominator = math.pow(phi_c2,2) * m * eta_LPT_star * eta_F_star
    x_opt = 1 / (1 + numerator / denominator)
    return x_opt


def calculate_P_spec(#норм
    q_T, v_take, m, phi_c1, phi_c2, x_opt, L_CB, eta_LPT_star, eta_F_star
):
    """
    Вычисляет параметр P_spec по формуле:

    P_yA = [(1 + q_T - v_take)/(m + 1)] * φ_c1 * √[2(1 - x_opt) * L_CB] +
            [m/(m + 1)] * φ_c2 * √[2(1 + q_T - v_take) * x_opt * L_CB * η*_THA * η*_B / m]

    Parameters:
    q_T (float): Параметр q_T
    v_take (float): Параметр v_take
    m (float): Параметр m
    phi_c1 (float): Коэффициент скорости φ_c1
    phi_c2 (float): Коэффициент скорости φ_c2
    x_opt (float): Оптимальный параметр x_opt
    L_CB (float): Параметр L_CB
    eta_LPT_star (float): КПД η*_THA
    eta_F_star (float): КПД η*_B

    Returns:
    float: Параметр P_yA
    """
    # Первое слагаемое
    term1 = ((1 + q_T - v_take) / (m + 1)) * phi_c1 * math.sqrt(2 * (1 - x_opt) * L_CB)

    # Второе слагаемое
    term2 = (
        (m / (m + 1))
        * phi_c2
        * math.sqrt(
            2 * (1 + q_T - v_take) * x_opt * L_CB * eta_LPT_star * eta_F_star / m
        )
    )

    P_yA = term1 + term2
    return P_yA

def calculate_l_free(#норм
    phi_c0,
    c_p,
    c_p_gas,
    t_g_star,
    pi_k_star,
    sigma_kc,
    sigma_bx,
    sigma_1,
    eta_p,
    t_h,
    q_t,
    nu_otb,
    eta_c,
    k,
    k_prime
):
    # Степени
    exponent_1 = (1 - k_prime) / k_prime
    exponent_2 = (k - 1) / k

    # Первая часть формулы
    bracket_inner = 1 - math.pow((pi_k_star * sigma_kc * sigma_bx * sigma_1), exponent_1)
    term1 = c_p_gas * t_g_star * bracket_inner * eta_p

    # Вторая часть формулы (вычитание)
    numerator = c_p * t_h * (math.pow(pi_k_star * sigma_bx, exponent_2) - 1)
    denominator = (1 + q_t - nu_otb) * eta_c
    term2 = numerator / denominator

    # Финальный результат
    result = (1 / math.pow(phi_c0, 2)) * (term1 - term2)

    return result

def calculate_c_spec(q_T, xi_intake, m, P_spec):#норм
    """
    Вычисляет удельный расход c_spec по формуле:
    
    c_spec = [3600 * q_T * (1 - ξ_orb)] / [(1 + m) * P_ȳ]
    
    Parameters:
    q_T (float): Параметр q_T
    xi_intake (float): Коэффициент потерь ξ_orb
    m (float): Параметр m
    P_spec (float): Параметр P_ȳ
    
    Returns:
    float: Удельный расход c_spec
    """
    c_c_spec = (3600 * q_T * (1 - xi_intake)) / ((1 + m) * P_spec)
    
    return c_c_spec

def calculate_phi_co(eta_T_star, pi_T_star, k_prime):#норм
    """
    Вычисляет коэффициент скорости φ_co по формуле:

    φ_co = 1 / [(1 - η*_T) * π*_T^((k'-1)/k') + η*_T]

    Parameters:
    eta_T_star (float): КПД η*_T
    pi_T_star (float): Степень понижения давления π*_T
    k_prime (float): Показатель адиабаты k'

    Returns:
    float: Коэффициент скорости φ_co
    """

    exponent = (k_prime - 1) / k_prime
    denominator = (1 - eta_T_star) * (math.pow(pi_T_star,exponent)) + eta_T_star
    phi_co = 1 / denominator

    return phi_co


def calculate_pi_c_star(k_prime):#норм
    """
    Вычисляет критическую степень повышения давления π_c* = π_cr по формуле:

    π_c* = [(k' + 1)/2]^(k'/(k'-1))

    Parameters:
    k_prime (float): Показатель адиабаты k'

    Returns:
    float: Критическая степень повышения давления π_c*
    """
    pi_c_star = math.pow(((k_prime + 1) / 2),(k_prime / (k_prime - 1)))

    return pi_c_star

def calculate_eta_c(sigma_bx, pi_k_star, k, eta_k_star):#норм
    exponent = (k - 1) / k

    # Числитель
    numerator = math.pow(sigma_bx * pi_k_star, exponent) - 1

    # Знаменатель
    term1 = math.pow(sigma_bx,exponent) * (math.pow(pi_k_star, exponent) - 1) * (1 / eta_k_star)
    term2 = math.pow(sigma_bx, exponent) - 1
    denominator = term1 + term2

    return numerator / denominator



def calculate_pi_t_star(sigma_bx, pi_k_star, sigma_kc, sigma_1, pi_c_star):#норм
    return (sigma_bx * pi_k_star * sigma_kc * sigma_1) / pi_c_star



def calculate_eta_p(pi_t_star, eta_t_star, pi_c_star, phi_c1, k_prime):#норм
    exponent = (1 - k_prime) / k_prime

    pi_t_pow = math.pow(pi_t_star, exponent)
    pi_c_pow = math.pow(pi_c_star, exponent)

    part1 = (1 - pi_t_pow) * eta_t_star
    part2 = (1 - part1) * (1 - pi_c_pow) * math.pow(phi_c1, 2)
    denominator = 1 - math.pow((pi_t_star * pi_c_star), exponent)

    return (part1 + part2) / denominator

def calc_u_cp(u_i,Dcp_i,D_i):
    return u_i * Dcp_i / D_i
def calc_sigma_p(material_density, ucp_i, d_rel_vt, f_rel):
    term1 = 2 * material_density * (ucp_i ** 2)
    
    term2 = (1 - d_rel_vt) / (1 + d_rel_vt)
    
    term3 = 1 - (1/3) * (1 - f_rel) * (1 + 1 / (1 + d_rel_vt))
    
    return term1 * term2 * term3
def select_stages(LTBD, effTBD, Dcp_tvd, Dcp_g, ucp_tvd,u_k1, D_vkvd):
    """
        z — подобранное число ступеней
        y — итоговый параметр Парсона
    """
    mu = 1.2
    c0 = math.sqrt(LTBD / effTBD * 2)   # теоретическая скорость истечения
    print(f"L - {LTBD:.4f}\nc0 - {c0:.4f}")
    z = 1                                         # начальное число ступеней
    deltaD = Dcp_tvd - Dcp_g
    stepD = 0# разница по диаметрам
    y = ucp_tvd / c0                                 # начальный параметр Парсона

    iter_count = 0
    opt = False

    while opt and iter_count < 10000:

        # корректируем число ступеней по диапазону
        if (y >= 0.6) and ((z - 1) > 0):
            z -= 1
        if y <= 0.45:
            z += 1
        stepD = deltaD / z

        # расчёт по ступеням
        for i in range(1, z + 1):
            u_cp_z =u_k1 * (Dcp_g + stepD * i) / D_vkvd
            y += math.pow(u_cp_z / c0, 2)

        # извлекаем корень
        y = math.sqrt(y)
        mu = LTBD / z / (u_cp_z * u_cp_z)
        # проверяем — подобрано?
        if (mu >= 1.2) and mu <= 1.8:
            opt = True 
        y = math.sqrt(effTBD/2/mu)
        if (y > 0.5) and (y < 0.6):
            opt = False
        else:
            opt = True
        
        iter_count += 1

    return {"z": z, "y": y, "mu": mu}
def calc_stages(LTBD, effTBD, Dcp_tvd, Dcp_g, ucp_tvd,u_k1, D_vkvd):
    c0 = math.sqrt(LTBD / effTBD * 2)   # теоретическая скорость истечения
    print(f"L - {LTBD:.4f}\nc0 - {c0:.4f}")
    z = 1                                         # начальное число ступеней
    y = ucp_tvd / c0                                 # начальный параметр Парсона
    mu = 0.0
    iter_count = 0
    opt = True
    while opt and iter_count < 10000:
        mu = LTBD / (z * (ucp_tvd ** 2)) 
        y = math.sqrt(effTBD / 2 /mu)
        if  1.2 <= mu <= 2.2 and 0.45 <= y <= 0.75:
            opt = False
        else:
            z+=1
        iter_count += 1

    return {"z": z, "y": y, "mu": mu}
def check_sigma(sigma_p, sigma_b):
        lower_bound = (1.5 + 2.2) / 2 * sigma_b
        res = sigma_p <= lower_bound
        if res:
            print(f"Допустимое значение σ_пр для ротора с σ_б={sigma_b} выполнено: σ_пр={sigma_p}")
        else:
            print(f"[ОШИБКА] Допустимое значение σ_пр для ротора с σ_б={sigma_b} НЕ выполнено: σ_пр={sigma_p}")
        return res