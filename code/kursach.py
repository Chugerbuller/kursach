import numpy as np
import math

def get_cp_air(T, R=287.0):
    poly = np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39])
    
    # Вычисление полинома (коэффициенты идут от старшей степени к младшей)
    cp = np.polyval(poly, T)
    
    cv = cp - R
    k = cp / cv
    
    return {'cp': cp, 'cv': cv, 'k': k}

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
            TK = TH * (1 + (pik**((kAir - 1) / kAir) - 1) / effK)
            
            # Сохраняем историю итераций
            history.append({
                'iteration': iter_count,
                'TK': TK,
                'cpAir': cpAir,
                'kAir': kAir
            })
            
            iter_count += 1
            
        except (ValueError, ZeroDivisionError) as e:
            print(f"Ошибка на итерации {iter_count}: {e}")
            break
    
    converged = abs(TK - TK_old) <= tol
    
    return {
        'TK': TK,
        'cpAir': cpAir,
        'cvAir': cvAir,
        'kAir': kAir,
        'iterations': iter_count,
        'converged': converged,
        'history': history
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
    R_CO2 = 8314.2 / (12 + 16 * 2)    # для углекислого газа
    R_H2O = 8314.2 / (2 + 16)         # для водяного пара
    R_N2 = 8314.2 / 28                # для азота
    R_O2 = 8314.2 / 32                # для кислорода
    
    # 2. Вычисляем низшую теплоту сгорания и теоретическое количество воздуха
    Hu = (33800 * C + 102500 * H) * 1000  # Дж/кг
    L0 = (8/3 * C + 8 * H) / 0.23         # кг/кг
    
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
            integral += coef * (T2**integrated_power - T1**integrated_power) / integrated_power
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
        g_CO2 = (11/3 * C) / (1 + alpha * L0)
        g_H2O = (9 * H) / (1 + alpha * L0)
        g_N2 = (0.77 * L0 * alpha) / (1 + L0 * alpha)
        g_O2 = (0.23 * (alpha - 1) * L0) / (1 + L0 * alpha)
        
        # Вычисляем средние теплоемкости в интервале [TK, TG]
        cp_CO2 = get_cp(TK, TG, np.array([-4.5e-7, 8.2e-4, -0.35, 850.0]))
        cp_H2O = get_cp(TK, TG, np.array([-3.8e-7, 7.1e-4, -0.32, 1860.0]))
        cp_N2 = get_cp(TK, TG, np.array([-2.8e-7, 6.5e-4, -0.28, 1040.0]))
        cp_O2 = get_cp(TK, TG, np.array([-3.1e-7, 6.8e-4, -0.29, 920.0]))
        #cp_air = get_cp(TK, TG, np.array([-3.2689e-7, 7.4230e-4, -3.1280e-1, 1042.39]))
        
        # Вычисляем теплоемкость и газовую постоянную смеси
        cp_mix = (cp_CO2 * g_CO2 + cp_H2O * g_H2O + 
                  cp_N2 * g_N2 + cp_O2 * g_O2)
        
        R_mix = (R_CO2 * g_CO2 + R_H2O * g_H2O + 
                 R_N2 * g_N2 + R_O2 * g_O2)
        
        cv_mix = cp_mix - R_mix
        
        # Обновляем показатель адиабаты
        k_old = k
        k = cp_mix / cv_mix
        
        # Обновляем коэффициент избытка воздуха
        alpha = 1/L0 * (Hu * effG / cp_mix / (TG - TK) - 1)
        
        # Сохраняем историю
        history.append({
            'iteration': iter_count,
            'alpha': alpha,
            'k': k,
            'cp_mix': cp_mix,
            'R_mix': R_mix,
            'g_CO2': g_CO2,
            'g_H2O': g_H2O,
            'g_N2': g_N2,
            'g_O2': g_O2
        })
        
        iter_count += 1
    
    # Проверяем сходимость
    converged = abs(k - k_old) <= tol
    
    return {
        'alpha': alpha,
        'k': k,
        'cp_mix': cp_mix,
        'R_mix': R_mix,
        'cv_mix': cv_mix,
        'mass_fractions': {
            'CO2': g_CO2,
            'H2O': g_H2O,
            'N2': g_N2,
            'O2': g_O2
        },
        'iterations': iter_count,
        'converged': converged,
        'history': history,
        'Hu': Hu,
        'L0': L0
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
        TB2 = TB * (1 + (piB**((kAir - 1) / kAir) - 1) / effB)
        
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
        cpGas, kGas, RGas, alpha = get_gas_parameters(tCorr, tCorr, Hu, L0, gC, effG, alpha)
        cpAir = get_cp_air(tCorr, tCorr)
        
        tCorrOld = tCorr
        
        cpMix = (cpGas * g1 + cpAir * g2) / (g1 + g2)
        RMix = (RGas * g1 + 287.0 * g2) / (g1 + g2)
        cvMix = cpMix - RMix
        kMix = cpMix / cvMix
        
        tCorr = (cpGasG * g1 * T1 + cpAirK * g2 * T2) / (cpGas * g1 + cpAir * g2)
        
        iter += 1
    
    return tCorr, cpMix, kMix, iter
def calculate_mixing_temperature(T1, T2, T3, g1, g2, g3, Hu, L0, gC, effG, alpha, R=287.0):
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
        cpGas, kGas, RGas, alpha = get_gas_parameters(tMix, tMix, Hu, L0, gC, effG, alpha)
        cpAir = get_cp_air(tMix, tMix)  # свойства воздуха при tMix
        
        tMixOld = tMix
        
        # Свойства смеси
        cpMix = (cpGas * g1 + cpAir * (g2 + g3)) / (g1 + g2 + g3)
        RMix = (RGas * g1 + R * (g2 + g3)) / (g1 + g2 + g3)
        cvMix = cpMix - RMix
        kMix = cpMix / cvMix
        
        # Новая температура смешения
        tMix = (cpGasG * g1 * T1 + cpAirK * g2 * T2 + cpAir2 * g3 * T3) / (cpGas * g1 + cpAir * (g2 + g3))
        
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
    
    switch = {
        0: case_0,
        1: case_1,
        2: case_2
    }
    
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
    deltaD = DCPTBD - DCPG  # различие в средних диаметрах между входным и выходным сечениями турбины
    stepD = 0  # инициализация переменной шага по диаметру
    y = uK1 / c0  # начальный параметр Парсонса (скорость на среднем диаметре / теоретическая скорость)
    
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

# Пример использования функции:
# z, y, iterations = calculate_turbine_stages(LTBD, effTBD, DCPTBD, DCPG, uK1, DBKBD)
# print(f"Оптимальное число ступеней: {z}, параметр Парсонса: {y:.3f}, итераций: {iterations}")