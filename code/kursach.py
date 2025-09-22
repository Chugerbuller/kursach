import numpy as np
import math


def get_cp_air(T, R=287.0):
    poly = np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39])

    # Вычисление полинома (коэффициенты идут от старшей степени к младшей)
    cp = np.polyval(poly, T)

    cv = cp - R
    k = cp / cv

    return {"cp": cp, "cv": cv, "k": k}


def get_cp_air(T1, T2):
    """
    Вычисление средней теплоемкости воздуха в интервале температур [T1, T2]
    через интегрирование аппроксимационного полинома

    Parameters:
    T1 : float
        Нижний предел температуры [K]
    T2 : float
        Верхний предел температуры [K]

    Returns:
    float : средняя изобарная теплоемкость в интервале [T1, T2] [Дж/(кг·К)]
    """
    # Коэффициенты интерполяционного полинома
    poly = np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39])
    n = len(poly)

    cp_integral = 0.0

    for i in range(n):
        # Степень полинома после интегрирования
        power = n - i  # оригинальная степень
        integrated_power = power + 1  # степень после интегрирования

        # Вычисление интеграла для каждого слагаемого
        term_T2 = poly[i] * T2**integrated_power / integrated_power
        term_T1 = poly[i] * T1**integrated_power / integrated_power
        cp_integral += term_T2 - term_T1

    # Средняя теплоемкость
    cp_avg = cp_integral / (T2 - T1)

    return cp_avg


def calculate_combustion_properties(
    C, H, L0, Hu, effG, TK, TG, cpCO2, cpH2O, cpN2, cpO2
):
    """
    Вычисляет свойства продуктов сгорания топлива

    Параметры:
    C - доля углерода в топливе
    H - доля водорода в топливе
    L0 - теоретическое количество воздуха
    Hu - низшая теплота сгорания топлива, Дж/кг
    effG - эффективность сгорания топлива
    TK - начальная температура, K
    TG - температура горения, K
    cpCO2, cpH2O, cpN2, cpO2 - удельные изобарные теплоемкости компонентов

    Возвращает:
    tuple: (k, alpha, gCO2, gH2O, gN2, gO2, R, cp, cv, iter_count)
    """

    # Вычисляем удельные газовые постоянные продуктов сгорания:
    RCO2 = 8314.2 / (12 + 16 * 2)  # для углекислого газа
    RH2O = 8314.2 / (2 + 16)  # для водяного пара
    RN2 = 8314.2 / 28  # для азота
    RO2 = 8314.2 / 32  # для кислорода

    # Задаем начальные значения
    k = 0
    alpha = 1
    kOld = 1.4
    iter_count = 0

    # Начинаем итерационный процесс
    while abs(k - kOld) > 1e-16 and iter_count < 100:
        # Вычисляем массовые доли компонентов
        gCO2 = 11 / 3 * C / (1 + alpha * L0)
        gH2O = 9 * H / (1 + alpha * L0)
        gN2 = 0.77 * L0 * alpha / (1 + L0 * alpha)
        gO2 = 0.23 * (alpha - 1) * L0 / (1 + L0 * alpha)

        # Вычисляем теплоемкости и газовую постоянную смеси
        cp = cpCO2 * gCO2 + cpH2O * gH2O + cpN2 * gN2 + cpO2 * gO2
        R = RCO2 * gCO2 + RH2O * gH2O + RN2 * gN2 + RO2 * gO2
        cv = cp - R

        kOld = k
        k = cp / cv

        # Вычисляем коэффициент избытка воздуха
        alpha = 1 / L0 * (Hu * effG / cp / (TG - TK) - 1)
        iter_count += 1

    return k, alpha, gCO2, gH2O, gN2, gO2, R, cp, cv, iter_count


def calculate_compressor_temperature(TH, pik, effK, R=287.0, max_iter=100, tol=1e-16):
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
            cpAir = get_cp_air(TH, TK)
            cvAir = cpAir - R

            # Защита от отрицательной теплоемкости
            if cvAir <= 0:
                raise ValueError("Отрицательная изохорная теплоемкость")

            kAir = cpAir / cvAir

            TK_old = TK
            TK = TH * (1 + (pik ** ((kAir - 1) / kAir) - 1) / effK)

            # Сохраняем историю итераций
            history.append(
                {"iteration": iter_count, "TK": TK, "cpAir": cpAir, "kAir": kAir}
            )

            iter_count += 1

        except (ValueError, ZeroDivisionError) as e:
            print(f"Ошибка на итерации {iter_count}: {e}")
            break

    converged = abs(TK - TK_old) <= tol

    return {
        "TK": TK,
        "cpAir": cpAir,
        "cvAir": cvAir,
        "kAir": kAir,
        "iterations": iter_count,
        "converged": converged,
        "history": history,
    }


def calculate_combustion_properties(C, H, TK, TG, effG, max_iter=1000, tol=1e-16):
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
    Hu = (33800 * C + 102500 * H) * 1000  # Дж/кг
    L0 = (8 / 3 * C + 8 * H) / 0.23  # кг/кг

    # 3. Функции для расчета теплоемкостей компонентов
    # (аналоги get_cp_air для каждого компонента)
    def get_cp(T1, T2, poly):
        """Теплоемкость CO2 в интервале [T1, T2]"""
        # Коэффициенты для CO2 (замените на реальные значения)
        n = len(poly)
        integral = 0.0
        for i, coef in enumerate(poly):
            power = n - i
            integrated_power = power + 1
            integral += (
                coef * (T2**integrated_power - T1**integrated_power) / integrated_power
            )
        return integral / (T2 - T1)

    # 4. Итерационный процесс
    k_old = 1.4  # начальное предположение
    k = 0.0
    alpha = 1.0  # начальный коэффициент избытка воздуха
    iter_count = 0
    converged = False

    history = []

    while abs(k - k_old) > tol and iter_count < max_iter:
        # Вычисляем массовые доли компонентов
        g_CO2 = (11 / 3 * C) / (1 + alpha * L0)
        g_H2O = (9 * H) / (1 + alpha * L0)
        g_N2 = (0.77 * L0 * alpha) / (1 + L0 * alpha)
        g_O2 = (0.23 * (alpha - 1) * L0) / (1 + L0 * alpha)

        # Вычисляем средние теплоемкости в интервале [TK, TG]
        cp_CO2 = get_cp(TK, TG, np.array([-4.5e-7, 8.2e-4, -0.35, 850.0]))
        cp_H2O = get_cp(TK, TG, np.array([-3.8e-7, 7.1e-4, -0.32, 1860.0]))
        cp_N2 = get_cp(TK, TG, np.array([-2.8e-7, 6.5e-4, -0.28, 1040.0]))
        cp_O2 = get_cp(TK, TG, np.array([-3.1e-7, 6.8e-4, -0.29, 920.0]))
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


def calculate_fan_temperature(TB, LB, effB, R=287.0):
    """Расчет температуры за вентилятором"""
    TB2old = 0.0
    TB2 = 300.0
    iter_count = 0

    while abs(TB2 - TB2old) > 1e-16 and iter_count < 100:
        cpAir = get_cp_air(TB, TB2)
        cvAir = cpAir - R
        kAir = cpAir / cvAir

        TB2old = TB2
        piB = (effB * LB / cpAir / TB + 1) ** (kAir / (kAir - 1))
        TB2 = TB * (1 + (piB ** ((kAir - 1) / kAir) - 1) / effB)

        iter_count += 1

    return TB2, piB, kAir, cpAir, iter_count


def get_gas_parameters(T, T_ref, Hu, L0, gC, effG, alpha):
    # Заглушка
    cpGas = 1200.0
    kGas = 1.3
    RGas = 280.0
    return cpGas, kGas, RGas, alpha


def calc_mix_temp(T1, T2, g1, g2, Hu, L0, gC, effG, alpha):
    # Истинные теплоемкости ДО цикла
    cpGasG = get_gas_parameters(T1, T1, Hu, L0, gC, effG, alpha)[0]
    cpAirK = get_cp_air(T2, T2)

    tCorr = T1
    tCorrOld = 0
    iter = 0

    while abs(tCorr - tCorrOld) > 1e-16 and iter < 100:
        cpGas, kGas, RGas, alpha = get_gas_parameters(
            tCorr, tCorr, Hu, L0, gC, effG, alpha
        )
        cpAir = get_cp_air(tCorr, tCorr)

        tCorrOld = tCorr

        cpMix = (cpGas * g1 + cpAir * g2) / (g1 + g2)
        RMix = (RGas * g1 + 287.0 * g2) / (g1 + g2)
        cvMix = cpMix - RMix
        kMix = cpMix / cvMix

        tCorr = (cpGasG * g1 * T1 + cpAirK * g2 * T2) / (cpGas * g1 + cpAir * g2)

        iter += 1

    return tCorr, cpMix, kMix, iter


def calculate_mixing_temperature(
    T1, T2, T3, g1, g2, g3, Hu, L0, gC, effG, alpha, R=287.0
):
    """
    Расчет температуры смешения трех потоков: газа, охлаждающего воздуха и воздуха второго контура
    """
    # Истинные теплоемкости ДО цикла (не меняются внутри цикла)
    cpGasG, kGas, RGas, alpha = get_gas_parameters(T1, T1, Hu, L0, gC, effG, alpha)
    cpAirK = get_cp_air(T2, T2)  # для охлаждающего воздуха при T2
    cpAir2 = get_cp_air(T3, T3)  # для воздуха второго контура при T3

    tMix = T1  # начальное приближение
    tMixOld = 0.0
    iter_count = 0

    while abs(tMix - tMixOld) > 1e-16 and iter_count < 100:
        # Свойства при текущей температуре смешения
        cpGas, kGas, RGas, alpha = get_gas_parameters(
            tMix, tMix, Hu, L0, gC, effG, alpha
        )
        cpAir = get_cp_air(tMix, tMix)  # свойства воздуха при tMix

        tMixOld = tMix

        # Свойства смеси
        cpMix = (cpGas * g1 + cpAir * (g2 + g3)) / (g1 + g2 + g3)
        RMix = (RGas * g1 + R * (g2 + g3)) / (g1 + g2 + g3)
        cvMix = cpMix - RMix
        kMix = cpMix / cvMix

        # Новая температура смешения
        tMix = (cpGasG * g1 * T1 + cpAirK * g2 * T2 + cpAir2 * g3 * T3) / (
            cpGas * g1 + cpAir * (g2 + g3)
        )

        iter_count += 1

    return tMix, cpMix, kMix, RMix, iter_count


def calculate_mix_properties(T1, T2, g1, g2, g3, Hu, L0, gC, effG, alpha, R=287.0):
    """
    Расчет свойств смеси газов без итерационного процесса
    """
    # Получаем средние теплоемкости в интервале [T1, T2]
    cpGas, kGas, RGas, alpha = get_gas_parameters(T1, T2, Hu, L0, gC, effG, alpha)
    cpAir = get_cp_air(T1, T2)

    # Свойства смеси по правилу смесей
    cpMix = (cpGas * g1 + cpAir * (g2 + g3)) / (g1 + g2 + g3)
    RMix = (RGas * g1 + R * (g2 + g3)) / (g1 + g2 + g3)
    cvMix = cpMix - RMix
    kMix = cpMix / cvMix

    return cpMix, RMix, cvMix, kMix


import math


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


def compressor_efficiency(pi_K_star, sigma_BX, k, eta_K_star):
    """
    Вычисляет КПД компрессора по формуле:

    η_c = [(σ_BX * π*_K)^((k-1)/k) - 1] /
          [σ_BX^((k-1)/k) * (π*_K^((k-1)/k) - 1) * (1/η*_K) + (σ_BX^((k-1)/k) - 1)]

    Parameters:
    pi_K_star (float): Степень повышения давления в компрессоре
    sigma_BX (float): Коэффициент потерь на входе
    k (float): Показатель адиабаты (отношение теплоемкостей cp/cv)
    eta_K_star (float): Изоэнтропический КПД компрессора

    Returns:
    float: КПД компрессора η_c
    """
    # Проверка входных параметров
    if pi_K_star <= 0:
        raise ValueError("Степень повышения давления должна быть > 0")
    if sigma_BX <= 0:
        raise ValueError("Коэффициент потерь на входе должен быть > 0")
    if k <= 1:
        raise ValueError("Показатель адиабаты должен быть > 1")
    if eta_K_star <= 0 or eta_K_star > 1:
        raise ValueError("Изоэнтропический КПД должен быть в диапазоне (0, 1]")

    # Вычисляем степени для упрощения формулы
    exponent = (k - 1) / k
    sigma_exp = sigma_BX**exponent
    pi_exp = pi_K_star**exponent

    # Числитель
    numerator = (sigma_BX * pi_K_star) ** exponent - 1

    # Знаменатель
    denominator = (sigma_exp * (pi_exp - 1) / eta_K_star) + (sigma_exp - 1)

    # КПД компрессора
    eta_c = numerator / denominator

    return eta_c


def optimal_parameter(phi_c1, phi_c2, q_T, v_take, m, eta_LPT_star, eta_F_star):
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
    numerator = phi_c1**2 * (1 + q_T - v_take)
    denominator = phi_c2**2 * m * eta_LPT_star * eta_F_star
    x_opt = 1 / (1 + numerator / denominator)
    return x_opt


def calculate_P_spec(
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
    import math

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


def calculate_full_free_energy(
    phi_c0,
    c_P_T_bar,
    T_G_star,
    pi_K_star,
    sigma_KC,
    sigma_BX,
    sigma_1,
    k_prime,
    eta_p,
    c_P_bar,
    T_H,
    k,
    q_T,
    v_take,
    eta_c,
):
    """
    Вычисляет параметр L_CB по формуле:

    L_CB = (1/φ_c0²) * [ c_P^T* T_F* (1 - (π_K* σ_KC σ_BX σ_1)^((1-k')/k')) η_p -
           (c_P T_H ((π_K* σ_BX)^((k-1)/k) - 1)) / ((1 + q_T - ν_or6) η_c) ]

    Parameters:
    phi_c0 (float): Коэффициент скорости φ_c0
    c_P_T_bar (float): Средняя удельная теплоемкость c_P^T*
    T_G_star (float): Температура T_F*
    pi_K_star (float): Степень повышения давления π_K*
    sigma_KC (float): Коэффициент σ_KC
    sigma_BX (float): Коэффициент σ_BX
    sigma_1 (float): Коэффициент σ_1
    k_prime (float): Показатель адиабаты k'
    eta_p (float): КПД η_p
    c_P_bar (float): Средняя удельная теплоемкость c_P
    T_H (float): Температура T_H
    k (float): Показатель адиабаты k
    q_T (float): Параметр q_T
    v_take (float): Параметр ν_or6
    eta_c (float): КПД η_c

    Returns:
    float: Параметр L_CB
    """

    # Первое слагаемое в скобках
    term1 = (
        c_P_T_bar
        * T_G_star
        * (1 - (pi_K_star * sigma_KC * sigma_BX * sigma_1) ** ((1 - k_prime) / k_prime))
        * eta_p
    )

    # Второе слагаемое в скобках
    term2 = (c_P_bar * T_H * ((pi_K_star * sigma_BX) ** ((k - 1) / k) - 1)) / (
        (1 + q_T - v_take) * eta_c
    )

    # Вычисляем L_CB
    L_CB = (1 / phi_c0**2) * (term1 - term2)

    return L_CB

def calculate_c_spec(q_T, xi_intake, m, P_spec):
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

def calculate_phi_co(eta_T_star, pi_T_star, k_prime):
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
    denominator = (1 - eta_T_star) * (pi_T_star**exponent) + eta_T_star
    phi_co = 1 / denominator

    return phi_co


def calculate_pi_c_star(k_prime):
    """
    Вычисляет критическую степень повышения давления π_c* = π_cr по формуле:

    π_c* = [(k' + 1)/2]^(k'/(k'-1))

    Parameters:
    k_prime (float): Показатель адиабаты k'

    Returns:
    float: Критическая степень повышения давления π_c*
    """
    pi_c_star = ((k_prime + 1) / 2) ** (k_prime / (k_prime - 1))

    return pi_c_star


def calculate_pi_T_star(sigma_BX, pi_K_star, sigma_KC, sigma_1, pi_c_star):
    """
    Вычисляет степень понижения давления π_T* по формуле:

    π_T* = (σ_BX * π_K* * σ_KC * σ_1) / π_c*

    Parameters:
    sigma_BX (float): Коэффициент σ_BX
    pi_K_star (float): Степень повышения давления π_K*
    sigma_KC (float): Коэффициент σ_KC
    sigma_1 (float): Коэффициент σ_1
    pi_c_star (float): Критическая степень повышения давления π_c*

    Returns:
    float: Степень понижения давления π_T*
    """
    pi_T_star = (sigma_BX * pi_K_star * sigma_KC * sigma_1) / pi_c_star

    return pi_T_star


def compression_efficiency(pi_K_star, sigma_BX, k, eta_K_star):
    """
    Вычисляет КПД сжатия по формуле:

    η_c = [(σ_BX * π*_K)^((k-1)/k) - 1] /
          [σ_BX^((k-1)/k) * (π*_K^((k-1)/k) - 1) * (1/η*_K) + (σ_BX^((k-1)/k) - 1)]

    Parameters:
    pi_K_star (float): Степень повышения давления в компрессоре
    sigma_BX (float): Коэффициент потерь на входе
    k (float): Показатель адиабаты (отношение теплоемкостей cp/cv)
    eta_K_star (float): Изоэнтропический КПД компрессора

    Returns:
    float: КПД компрессора η_c
    """
    # Проверка входных параметров
    if pi_K_star <= 0:
        raise ValueError("Степень повышения давления должна быть > 0")
    if sigma_BX <= 0:
        raise ValueError("Коэффициент потерь на входе должен быть > 0")
    if k <= 1:
        raise ValueError("Показатель адиабаты должен быть > 1")
    if eta_K_star <= 0 or eta_K_star > 1:
        raise ValueError("Изоэнтропический КПД должен быть в диапазоне (0, 1]")

    # Вычисляем степени для упрощения формулы
    exponent = (k - 1) / k
    sigma_exp = sigma_BX**exponent
    pi_exp = pi_K_star**exponent

    # Числитель
    numerator = (sigma_BX * pi_K_star) ** exponent - 1

    # Знаменатель
    denominator = (sigma_exp * (pi_exp - 1) / eta_K_star) + (sigma_exp - 1)

    # КПД сжатия
    eta_c = numerator / denominator

    return eta_c


def expansion_efficiency(
    pi_K_star, sigma_BX, sigma_KC, k_prime, eta_T_star, sigma_1, phi_c1
):
    """
    Вычисляет КПД расширения по формуле:

    η_p = [(1-π*_T)^((1-k')/k') * η*_T + (1-(1-π*_T)^((1-k')/k')) * η*_T) * (1-π*_C)^((1-k')/k') * phi_c1^2]/
          [1-(π*_T * π*_C)^((1-k')/k')]

    Parameters:
    pi_K_star (float): Степень повышения давления в компрессоре
    sigma_BX (float): Коэффициент потерь на входе
    sigma_KC (float): Коэффициент потерь в камере сгорания
    k_prime (float): Показатель адиабаты k'
    eta_T_star (float): Изоэнтропический КПД турбины

    Returns:
    float: КПД расширения η_p
    """
    # Проверка входных параметров
    if pi_K_star <= 0:
        raise ValueError("Степень повышения давления должна быть > 0")
    if sigma_BX <= 0:
        raise ValueError("Коэффициент потерь на входе должен быть > 0")
    if sigma_KC <= 0:
        raise ValueError("Коэффициент потерь в камере сгорания должен быть > 0")
    if k_prime <= 1:
        raise ValueError("Показатель адиабаты должен быть > 1")
    if eta_T_star <= 0 or eta_T_star > 1:
        raise ValueError("Изоэнтропический КПД должен быть в диапазоне (0, 1]")

    # Вычисляем степени для упрощения формулы
    exponent = (1 - k_prime) / k_prime
    pi_T_exp = pi_T_star**exponent
    pi_C_exp = pi_C_star**exponent
    pi_exp = (pi_T_star * pi_C_exp) ** exponent

    # Вычисляем степень расширения в сопле. Режим расчётный
    pi_C_star = ((k_prime + 1) / 2) ** exponent

    # Вычисляем степень понижения давления в турбине
    pi_T_star = sigma_BX * pi_K_star * sigma_KC * sigma_1 / pi_C_star

    numerator = (1 - pi_T_exp) * eta_T_star + (1 - (1 - pi_T_exp) * eta_T_star) * (
        1 - pi_C_exp
    ) * phi_c1**2
    denominator = 1 - pi_exp
    eta_p = numerator / denominator
    return eta_p


# Пример использования функции:
# z, y, iterations = calculate_turbine_stages(LTBD, effTBD, DCPTBD, DCPG, uK1, DBKBD)
# print(f"Оптимальное число ступеней: {z}, параметр Парсонса: {y:.3f}, итераций: {iterations}")
