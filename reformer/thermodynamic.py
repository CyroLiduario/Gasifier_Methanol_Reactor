import numpy as np

def capacidade_calorifica_gases(componentes_gases, Tg):
    tau = max(0.0, Tg) / 1000.0
    cp = {}
    for comp, dados in componentes_gases.items():
        c = dados["cp_coefs"]
        # O (MM/1000) de conversão para kg/mol corta com o 1000 da fórmula original
        cp[comp] = dados["MM"] * (c[0] + c[1]*tau + c[2]*(tau**2) + c[3]*(tau**3))
    return cp

def derivada_capacidade_calorifica_gases(componentes_gases, Tg):
    tau = max(0.0, Tg) / 1000.0
    dcp_dT = {}
    for comp, dados in componentes_gases.items():
        c = dados["cp_coefs"]
        # Ajuste exato da derivada analítica do Cyro
        dcp_dT[comp] = (dados["MM"] / 1000.0) * (c[1] + 2.0*c[2]*tau + 3.0*c[3]*(tau**2))
    return dcp_dT

def capacidade_calorifica_solidos(componentes_solidos, Ts):
    cp_s = {}
    # Biomassa e Umidade: o fator de 1000 da fórmula corta com o g/mol do JSON
    dados_b = componentes_solidos["Biomassa"]
    cp_s["Biomassa"] = dados_b["MM"] * dados_b["cp_coefs"][0]
    
    dados_u = componentes_solidos["Umidade"]
    cp_s["Umidade"] = dados_u["MM"] * dados_u["cp_coefs"][0]
    
    # Char: A fórmula original do Cyro *NÃO* tinha o multiplicador 10^3
    dados_c = componentes_solidos["Char"]
    c_c = dados_c["cp_coefs"]
    t_celsius = max(0.0, Ts) - 273.15
    cp_s["Char"] = (dados_c["MM"] / 1000.0) * (c_c[0] + c_c[1]*t_celsius + c_c[2]*(t_celsius**2))
    
    return cp_s

def derivada_capacidade_calorifica_solidos(componentes_solidos, Ts):
    dcp_s_dT = {"Biomassa": 0.0, "Umidade": 0.0}
    dados_c = componentes_solidos["Char"]
    c_c = dados_c["cp_coefs"]
    t_celsius = max(0.0, Ts) - 273.15
    dcp_s_dT["Char"] = (dados_c["MM"] / 1000.0) * (c_c[1] + 2.0*c_c[2]*t_celsius)
    return dcp_s_dT

def viscosidade_gas(Tg):
    """
    Calcula a viscosidade dinâmica da mistura gasosa (Pa.s) simplificada pelo modelo.
    """
    return 1.98e-5 * ((max(1e-10, Tg) / 300.0)**(2.0 / 3.0))

def propriedades_transporte(Tg, Ts, rhom_g, U_g, mu, cpm_g, reator_params):
    """
    Calcula todos os parâmetros de transferência de calor convectivos e radiativos
    entre as fases gás, sólido e a parede do reator.
    
    Retorna um dicionário com: Qsg, Qsw, Qgw (Fluxos de calor) e coeficientes auxiliares.
    """
    # Parâmetros estruturais e operacionais do reator
    D = reator_params["D"]
    eps = reator_params["eps"]
    phi = reator_params["phi"]
    em = reator_params["em"]
    d_p = reator_params["d_p"]
    T_w = reator_params["T_w"]
    sigma = reator_params["sigma"]
    A_p = 6.0 * (1.0 - eps) / d_p
    
    # Evitar divisões por zero ou valores físicos impossíveis nas raízes
    Tg_safe = max(1e-10, Tg)
    Ts_safe = max(1e-10, Ts)
    U_g_safe = max(1e-10, U_g)
    
    # 1. Condutividades térmicas puras
    lambda_g = 4.8e-4 * (Tg_safe**0.717)
    lambda_s = 0.13 + 0.0003 * (Ts_safe - 273.15)
    kappa = lambda_s / lambda_g
    
    # 2. Números adimensionais
    Re = max(1e-10, d_p * U_g_safe * rhom_g / mu)
    Pr = max(1e-10, cpm_g * mu / lambda_g)
    
    # 3. Coeficientes de transferência por radiação
    hrv = 4.0 * sigma / (1.0 + 0.5 * (eps / (1.0 - eps)) * ((1.0 - em) / em)) * (Tg_safe**3)
    hrs = 4.0 * sigma * (em / (2.0 - em)) * (Ts_safe**3)
    
    # 4. Condutividades térmicas efetivas do leito (Modelo de Hobbs)
    lambdarad0 = (eps * lambda_g * (1.0 + d_p * hrv / lambda_g) + 
                  lambda_g * (1.0 - eps) / (((1.0 / phi + d_p * hrs / lambda_g)**(-1)) + 2.0 / (3.0 * kappa)))
    
    lambdarad_g = lambda_g * (eps * (1.0 + d_p * hrv / lambda_g) + 0.14 * Re * Pr / (1.0 + 46.0 * ((d_p / D)**2)))
    lambdarad_s = lambda_g * (1.0 - eps) * (((1.0 / phi + d_p * hrs / lambda_g)**(-1)) + 2.0 / (3.0 * kappa))**(-1)
    
    # 5. Coeficientes globais de transferência de calor na parede
    hw = (2.44 * lambdarad0 / (D**(4.0 / 3.0))) + (0.033 * lambda_g / d_p) * Re * Pr
    hsw = (hw * lambdarad_s) / (lambdarad_g + lambdarad_s)
    hgw = (hw * lambdarad_g) / (lambdarad_g + lambdarad_s)
    
    # 6. Coeficiente de transferência de calor Sólido-Gás (Convecção Interfásica)
    hsg = (2.06 * cpm_g * rhom_g * U_g_safe / eps) * (Re**(-0.575)) * (Pr**(-2.0 / 3.0))
    
    # 7. Fluxos volumétricos/lineares de calor (Termos fonte das EDOs de energia)
    Qsg = hsg * A_p * (Ts - Tg)
    Qsw = (hsw * 10.0 / D) * (Ts - T_w)
    Qgw = (hgw * 10.0 / D) * (Tg - T_w)
    
    return {
        "Qsg": Qsg,
        "Qsw": Qsw,
        "Qgw": Qgw,
        "hsg": hsg,
        "hw": hw,
        "Re": Re,
        "Pr": Pr
    }