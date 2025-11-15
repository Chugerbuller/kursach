import kursach as kurs
import matplotlib.pyplot as plt


class Table_temperature:
    def __init__(self, t, m):
        super().__init__()
        self.t = t
        self.m = m

class Table_m:
    def __init__(self, m, pi_k_full):
        super().__init__()
        self.m = m
        self.pi_k_full = pi_k_full


class Table_pi_k_full:
    def __init__(
        self,
        pi_k_full,
        t_gas_full,
        alpha,
        x_opt,
        k,
        k_gas,
        p_spec,
        c_spec,
        l_free_energy,
        eff_comp,
        T_k,
    ):
        super().__init__()
        self.pi_k_full = pi_k_full
        self.t_gas_full = t_gas_full
        self.alpha = alpha
        self.x_opt = x_opt
        self.k = k
        self.k_gas = k_gas
        self.p_spec = p_spec
        self.c_spec = c_spec
        self.l_free_energy = l_free_energy
        self.eff_comp = eff_comp
        self.T_k = T_k

def calc_opt_params(engine, coef):
    t_array = []
    m_array = []
    pi_array = []
    Pik_full = []
    
    T_gas_full = [
        engine["T_gas_full"] - 150,
        engine["T_gas_full"],
        engine["T_gas_full"] + 150,
    ]
    m = [engine["m"] * 0.8,
        engine["m"],
        engine["m"] * 1.2]
    
    step = 0.8
    for i in range(7):
        Pik_full.append(step * engine["Pik_full"])
        step += 0.05

    for T_gas_full_i in T_gas_full:
        for m_i in m:
            for Pik_full_i in Pik_full:
                pi_array.append(calc_proto(coef,T_gas_full_i,m_i,Pik_full_i))
            m_array.append(Table_m(m_i, pi_array.copy()))
            pi_array.clear()
        t_array.append(Table_temperature(T_gas_full_i, m_array.copy()))
        m_array.clear()

    return t_array

def calc_proto(coef, T_gas_full_i, m_i, Pik_full_i):
    
    compressor = kurs.calculate_compressor_temperature(
                288, Pik_full_i, coef["effk_comp_full"]
            )
    T_k = compressor['TK']
    k_air = compressor['kAir']
    cp_air = compressor['cpAir']
    combustion_props = kurs.calculate_combustion_properties(
                0.86, 0.14, T_k, T_gas_full_i, coef["effk_gas"]
            )
    cp_bar = combustion_props["cp_mix"]
    k_gas = combustion_props["k"]
    q_T = 1 / (combustion_props["alpha"] * combustion_props["L0"])
    v_take = coef["ksi_take"] - coef["g_air_back"]
    phi_c_star = kurs.calculate_pi_c_star(k_gas)
    comp_eff = kurs.calculate_eta_c(coef['sigma_intake'],
                                    Pik_full_i,
                                    k_air,
                                    coef['effk_comp_full'])
    phi_t_star = kurs.calculate_pi_t_star(coef['sigma_intake'],
                                            Pik_full_i,
                                            coef['sigma_cc'],
                                            coef['sigma_1'],
                                            phi_c_star)
    eta_t_star = ((coef['effk_hpt_full'] + coef['effk_lpt_full']) / 2) * (1 + 0.03)
    expansion_eff = kurs.calculate_eta_p(phi_t_star,eta_t_star,phi_c_star,coef['phi_c1'],k_gas)
    phi_c0 = kurs.calculate_phi_co(eta_t_star,phi_t_star,k_gas)
    L_free = kurs.calculate_l_free(
        phi_c0,
        cp_air,
        cp_bar,
        T_gas_full_i,
        Pik_full_i,
        coef['sigma_cc'],
        coef['sigma_intake'],
        coef['sigma_1'],
        expansion_eff,
        288,
        q_T,
        v_take,
        comp_eff,
        k_air,
        k_gas
        )
    x_opt = kurs.optimal_parameter(
        coef['phi_c1'],
        coef['phi_c2'],
        q_T,
        v_take,
        m_i,
        coef['effk_lpt_full'],
        coef['effk_fan_full']
        )
    p_spec = kurs.calculate_P_spec(
                    q_T,
                    v_take,
                    m_i,
                    coef["phi_c1"],
                    coef["phi_c2"],
                    x_opt,
                    L_free,
                    coef["effk_lpt_full"],
                    coef["effk_fan_full"],
                )
    c_spec = kurs.calculate_c_spec(q_T, coef["ksi_take"], m_i, p_spec)
    return Table_pi_k_full(
                        Pik_full_i,
                        T_gas_full_i,
                        combustion_props["alpha"],
                        x_opt,
                        k_air,
                        k_gas,
                        p_spec,
                        c_spec,
                        L_free,
                        comp_eff,
                        T_k
                    )
                