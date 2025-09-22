import kursach as kurs


def calc_opt_params():

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

    T_gas_full = [
        engine["T_gas_full"] - 150,
        engine["T_gas_full"],
        engine["T_gas_full"] + 150,
    ]
    m = [engine["m"] * 0.8, engine["m"], engine["m"] * 1.2]
    Pik_full = []
    step = 0.8
    for i in range(7):
        Pik_full.append(step * engine["Pik_full"])
        step += 0.05
    for T_gas_full_i in T_gas_full:
        for m_i in m:
            for Pik_full_i in Pik_full:
                T_k = kurs.calculate_compressor_temperature(
                    T_gas_full_i, Pik_full_i, coef["effk_comp_full"]
                )["TK"]
                combustion_props = kurs.calculate_combustion_properties(
                    0.86, 0.14, T_k, T_gas_full_i, coef["effk_gas"]
                )
                k_gas = combustion_props["k"]
                q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
                v_take = coef["ksi_take"] - coef["g_air_back"]
                x_opt = kurs.optimal_parameter(
                    coef["phi_c1"],
                    coef["phi_c2"],
                    q_T,
                    v_take,
                    m_i,
                    coef["effk_lpt_full"],
                    coef["effk_fan_full"],
                )
                phi_c_star = kurs.calculate_pi_c_star(k_gas)
                phi_t_star = kurs.calculate_pi_T_star(
                    coef["sigma_intake"],
                    Pik_full_i,
                    coef["sigma_cc"],
                    coef["sigma_1"],
                    phi_c_star,
                )
                L_free = kurs.calculate_full_free_energy(
                    kurs.calculate_phi_co(
                        (coef["effk_hpt_full"] * coef["effk_lpt_full"]),
                        phi_t_star,
                        k_gas,
                    ),
                    combustion_props["cp_mix"],
                    T_gas_full_i,
                    Pik_full_i,
                    coef["sigma_cc"],
                    coef["sigma_intake"],
                    coef["sigma_1"],
                    k_gas,
                    kurs.expansion_efficiency(
                        Pik_full_i,
                        coef["sigma_intake"],
                        coef["sigma_cc"],
                        k_gas,
                        coef["effk_hpt_full"] * coef["effk_lpt_full"],
                        coef["sigma_1"],
                        coef["phi_c1"],
                    ),
                    kurs.get_cp_air(T_gas_full_i, 288),
                    288,
                    1.4,
                    q_T,
                    v_take,
                    kurs.compressor_efficiency(
                        Pik_full_i, coef["sigma_intake"], 1.4, coef["effk_comp_full"]
                    ),
                )
