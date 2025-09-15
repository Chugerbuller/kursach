import kursach as kurs
 
def calc_opt_params()

    engine = {
        'P': 264.447
        'Pik_full': 34.5
        'G_air': 727.56
        'm': 4.10
        'T_gas_full': 1631
        'C_spec': 0.0357  
    }

    coef = {
        'g_c': 0.86
        'sigma_intake': 0.99
        'sigma_cc': 0.955
        'sigma_1': 0.99
        'phi_c1': 0.98
        'phi_c2': 0.98
        'eta_gas': 0.995
        'eta_comp_full': 0.81
        'eta_fan_full': 0.87
        'a': 0.03
        'eta_hpt_full': 0.88
        'eta_lpt_full': 0.91
        'ksi_take': 0.155
        'g_air_back': 0.124
    }

    T_gas_full = [engine['T_gas_full'] - 150, engine['T_gas_full'], engine['T_gas_full'] + 150]
    Pik_full = []
    step = 0.8
    for i in range(7):
        Pik_full.append(step*engine['Pik_full'])
        step += 0.05
    m = [engine['m'] * 0.8, engine['m'], engine['m'] * 1.2]
    