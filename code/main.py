import lab1 as lab
import matplotlib.pyplot as plt
import numpy as np # Эта библиотека для работы с числами и массивами
import json

engine = {
        "P": 264.447,
        "Pik_full": 34.5,
        "G_air": 727.56,
        "m": 4.10,
        "T_gas_full": 1631,
        "C_spec": 0.0357,
    }
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

X = 1.07
p_project = engine["P"] * X
proto = lab.calc_proto(coef,engine["T_gas_full"], engine["m"], engine["Pik_full"])
l_free_project = proto.l_free_energy * X * X

Pik_full = []
T_gas_full = [
        engine["T_gas_full"] - 150,
        engine["T_gas_full"],
        engine["T_gas_full"] + 150,
        engine["T_gas_full"] + 200,
    ]
m = [engine["m"] * 0.8,
        engine["m"],
         engine["m"] * 1.2]
    
step = 0.4
t_temp = []
m_temp = []

for i in range(18):
        Pik_full.append(step * engine["Pik_full"])
        step += 0.05
for T_gas_full_i in T_gas_full:
    
    m_temp.append(lab.Table_m(engine["m"],
                                lab.calc_proto(coef
                                                ,T_gas_full_i
                                                ,engine["m"],
                                                engine["Pik_full"])))
    t_temp.append(lab.Table_temperature(T_gas_full_i, m_temp.copy()))
    m_temp.clear()
min_energy = 0
max_energy = 9999999999999
min_t = 0
max_t = 0
for t_table in t_temp:
    free_energy =  t_table.m[0].pi_k_full.l_free_energy
    if free_energy > min_energy and free_energy < l_free_project:
        min_energy = free_energy
        min_t = t_table.t
    if free_energy < max_energy and free_energy > l_free_project:
        max_energy = free_energy
        max_t = t_table.t
T_gas_opt = min_t + (((max_t - min_t) / (max_energy - min_energy)) * (l_free_project - min_energy))

opt_pi_map = {}

for pi_k_full_i in Pik_full:
    opt_pi_map[pi_k_full_i] = lab.calc_proto(coef,T_gas_opt, engine["m"], pi_k_full_i)
max_eff = 0.0
opt_pi = 0.0
for pi in opt_pi_map:
    if opt_pi_map[pi].eff_comp > max_eff:
        max_eff = opt_pi_map[pi].eff_comp
        opt_pi = pi

print(opt_pi)
fuel_map = {}
min_fuel = 99999999
opt_m = 0.0

for m_i in m:
    fuel_map[m_i] = lab.calc_proto(coef,T_gas_opt, m_i, opt_pi)
    
for m_i in fuel_map:
    if fuel_map[m_i].c_spec < min_fuel:
        min_fuel = fuel_map[m_i].c_spec
        opt_m = m_i
        
# Создадим данные для примера
# Создаем массив чисел от 0 до 10 с шагом 0.1
y = []
test = []
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, opt_m, pi_i).c_spec)
    test.append(lab.calc_proto(coef,T_gas_opt, opt_m, pi_i))
temp = lab.calc_proto(coef,engine["T_gas_full"], engine['m'], engine["Pik_full"])
# Создаем график
plt.figure(figsize=(15, 10)) # Задаем размер картинки (ширина, высота)
plt.plot(Pik_full, y, label='opt', color='blue', linewidth=1) # Рисуем линию с меткой и цветом
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,engine["T_gas_full"], engine['m'], pi_i).c_spec)
plt.plot(Pik_full, y, label='proto', color='black', linewidth=1)
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, 4.8, pi_i).c_spec)
plt.plot(Pik_full, y, label='var', color='red', linewidth=1)
# Добавляем названия и легенду
plt.title('График удельного расхода топлива') # Заголовок
plt.xlabel('πK*') # Подпись оси X
plt.ylabel('Удельный расход топлива, кг/Н/ч') # Подпись оси Y
plt.grid(True) # Включаем сетку
plt.legend() # Показываем легенду

# Самое важное: сохраняем график в файл
# dpi (dots per inch) отвечает за качество картинки
plt.savefig('c-spec.png', dpi=300, bbox_inches='tight')

y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, opt_m, pi_i).p_spec)

# Создаем график
plt.figure(figsize=(15, 10)) # Задаем размер картинки (ширина, высота)
plt.plot(Pik_full, y, label='opt', color='blue', linewidth=1) # Рисуем линию с меткой и цветом
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,engine["T_gas_full"], engine['m'], pi_i).p_spec)
plt.plot(Pik_full, y, label='proto', color='black', linewidth=1)
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, 4.8, pi_i).p_spec)
plt.plot(Pik_full, y, label='var', color='red', linewidth=1)

# Добавляем названия и легенду
plt.title('График удельной тяги') # Заголовок
plt.xlabel('πK*') # Подпись оси X
plt.ylabel('Удельная тяга, м/с ') # Подпись оси Y
plt.grid(True) # Включаем сетку
plt.legend() # Показываем легенду

# Самое важное: сохраняем график в файл
# dpi (dots per inch) отвечает за качество картинки
plt.savefig('p-spec.png', dpi=300, bbox_inches='tight')
print()
res = lab.calc_opt_params(engine, coef)

y.clear()
plt.figure(figsize=(15, 10))

for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, opt_m, pi_i).l_free_energy)
plt.plot(Pik_full, y, label='opt', color='blue', linewidth=1) # Рисуем линию с меткой и цветом # Рисуем линию с меткой и цветом
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,engine["T_gas_full"], engine['m'], pi_i).l_free_energy)
plt.plot(Pik_full, y, label='proto', color='black', linewidth=1)
y.clear()
for pi_i in Pik_full:
    y.append(lab.calc_proto(coef,T_gas_opt, 4.8, pi_i).l_free_energy)
plt.plot(Pik_full, y, label='var', color='red', linewidth=1)
# Добавляем названия и легенду
plt.title('График свободной энергии') # Заголовок
plt.xlabel('πK*') # Подпись оси X
plt.ylabel('свободная энергия дж/кг') # Подпись оси Y
plt.grid(True) # Включаем сетку
plt.legend() # Показываем легенду

# Самое важное: сохраняем график в файл
# dpi (dots per inch) отвечает за качество картинки
plt.savefig('l-free.png', dpi=300, bbox_inches='tight')

res = lab.calc_opt_params(engine,coef)
res_string = "T\tm\t\tpi|T*|a|Xopt|k|k_gas|P_cpec|C_cpec|L_free\n"
for T_table in res:
    res_string += (f"T:{T_table.t:.1f}\n")
    for m_table in T_table.m:
        res_string += (f"\tm:{m_table.m:.2f}\n")
        for pi_table in m_table.pi_k_full:
            res_string += (f"\t\t{pi_table.pi_k_full:.4f}|{pi_table.t_gas_full:.4f}|{pi_table.alpha:.4f}|{pi_table.x_opt:.4f}|{pi_table.k:.4f}|{pi_table.k_gas:.4f}|{pi_table.p_spec:.4f}|{pi_table.c_spec:.4f}|{pi_table.l_free_energy:.4f}\n")
print(res_string)
with open("results.txt", "w") as file:  
    file.write(res_string)  