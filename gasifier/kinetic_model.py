import math
from thermodynamic import (
    capacidade_calorifica_gases, derivada_capacidade_calorifica_gases,
    capacidade_calorifica_solidos, derivada_capacidade_calorifica_solidos,
    viscosidade_gas, propriedades_transporte
)

def gasifier_odes(z, y, materiais, reator_params):
    """
        Calcula as derivadas espaciais das variáveis de estado
        e as taxas de reação do gaseificador.

        Retorna:
            dydz: vetor com as 16 derivadas espaciais.
            taxas_reacao: dicionário com as taxas instantâneas das reações.
    """
    # 1. Desempacotamento das variáveis de estado (y)
    rho_B, rho_char, rho_M = y[0], y[1], y[2]
    C_g_H2O, C_g_O2, C_g_CO, C_g_CO2, C_g_H2, C_g_CH4, C_g_N2, C_g_Tar = y[3:11]
    U_g, P, U_s, T_s, T_g = y[11:16]

    # Prevenção numérica contra valores negativos que quebram logaritmos/raízes
    T_s = max(1e-4, T_s)
    T_g = max(1e-4, T_g)
    
    # Constante universal dos gases
    Rg = 8.314 # J/(mol.K)
    CRF = 1000.0 # Fator de ajuste cinético original do modelo
    
    # 2. Concentrações, Densidades e Frações
    gases_keys = ["H2O", "O2", "CO", "CO2", "H2", "CH4", "N2", "Tar"]
    C_g = [max(0, c) for c in y[3:11]]
    Cm_g = sum(C_g) + 1e-10 # Concentração molar total
    
    rho_g = {}
    for i, comp in enumerate(gases_keys):
        rho_g[comp] = C_g[i] * materiais["Gases"][comp]["MM"] / 1000.0 # kg/m³
        
    rhom_g = sum(rho_g.values())
    rhom_s = max(0, rho_B) + max(0, rho_char) + max(0, rho_M)
    
    # Fator de existência de carvão (char) para ativar reações gás-sólido
    fatorchar = 1.0 if rho_char > 1e-4 else 0.0

    # 3. Termodinâmica e Transferência de Calor
    cp_g = capacidade_calorifica_gases(materiais["Gases"], T_g)
    dcp_g_dT = derivada_capacidade_calorifica_gases(materiais["Gases"], T_g)
    
    cp_s = capacidade_calorifica_solidos(materiais["Solidos"], T_s)
    dcps_dT = derivada_capacidade_calorifica_solidos(materiais["Solidos"], T_s)
    
    mu = viscosidade_gas(T_g)
    
    # Cp médio da mistura gasosa
    cpm_g = sum((cp_g[comp] / (materiais["Gases"][comp]["MM"]/1000.0)) * rho_g[comp] for comp in gases_keys) / rhom_g
    
    transporte = propriedades_transporte(T_g, T_s, rhom_g, U_g, mu, cpm_g, reator_params)
    
    # 4. Pressões Parciais e Constantes de Equilíbrio (K)
    p_parcial = {} # Pressões parciais em bar
    for i, comp in enumerate(gases_keys):
        p_parcial[comp] = max(0, P * (C_g[i] / Cm_g) / 1e5)
        
    K_eq = materiais["Cinetica"]["Constantes_Equilibrio"]
    K1 = 10**((K_eq["K1"]["A"] / T_s) + K_eq["K1"]["B"])
    K2 = 10**((K_eq["K2"]["A"] / T_s) + K_eq["K2"]["B"])
    K3 = 10**((K_eq["K3"]["A"] / T_s) + K_eq["K3"]["B"])
    K4 = 10**((K_eq["K4"]["A"] / T_g) + K_eq["K4"]["B"])
    K5 = 10**((K_eq["K5"]["A"] / T_g) + K_eq["K5"]["B"])
    K6 = 10**((K_eq["K6"]["A"] / T_g) + K_eq["K6"]["B"])
    K7 = 10**((K_eq["K7"]["A"] / T_g) + K_eq["K7"]["B"])

    # 5. Taxas de Reação Cinética (Arrhenius)
    cinetica = materiais["Cinetica"]["Reacoes"]
    
    # Secagem e Pirólise
    m = max(0, rho_M) * cinetica["secagem_m"]["A"] * math.exp(-cinetica["secagem_m"]["E"] / (Rg * T_s))
    Rp1 = max(0, rho_B) * cinetica["pirolise_p1"]["A"] * math.exp(-cinetica["pirolise_p1"]["E"] / (Rg * T_s))
    
    # Reações heterogêneas (Gás-Sólido) - Dependentes de T_s
    R1 = fatorchar * CRF * cinetica["g1_boudouard"]["A"] * math.exp(-cinetica["g1_boudouard"]["E"] / (Rg * T_s)) * (p_parcial["CO2"] - (p_parcial["CO"]**2) / K1)
    R2 = fatorchar * CRF * cinetica["g2_gasificacao_H2O"]["A"] * math.exp(-cinetica["g2_gasificacao_H2O"]["E"] / (Rg * T_s)) * (p_parcial["H2O"] - (p_parcial["CO"] * p_parcial["H2"]) / K2)
    R3 = fatorchar * CRF * cinetica["g3_gasificacao_H2"]["A"] * math.exp(-cinetica["g3_gasificacao_H2"]["E"] / (Rg * T_s)) * ((p_parcial["H2"]**2) - p_parcial["CH4"] / K3)

    # Reações homogêneas (Gás-Gás) - Dependentes de T_g
    R4 = CRF * cinetica["g4_reformacao_CH4"]["A"] * math.exp(-cinetica["g4_reformacao_CH4"]["E"] / (Rg * T_g)) * ((p_parcial["CH4"] * p_parcial["H2O"]) - (p_parcial["CO"] * (p_parcial["H2"]**3)) / K4)
    R5 = CRF * cinetica["c1_oxidacao_CH4"]["A"] * math.exp(-cinetica["c1_oxidacao_CH4"]["E"] / (Rg * T_g)) * ((p_parcial["CH4"] * (p_parcial["O2"]**2)) - (p_parcial["CO2"] * (p_parcial["H2O"]**2)) / K5)
    R6 = CRF * cinetica["c2_oxidacao_CO"]["A"] * math.exp(-cinetica["c2_oxidacao_CO"]["E"] / (Rg * T_g)) * ((p_parcial["CO"]**2) * p_parcial["O2"] - (p_parcial["CO2"]**2) / K6)
    R7 = CRF * cinetica["c3_oxidacao_H2"]["A"] * math.exp(-cinetica["c3_oxidacao_H2"]["E"] / (Rg * T_g)) * ((p_parcial["H2"]**2) * p_parcial["O2"] - (p_parcial["H2O"]**2) / K7)

    # Water-Gas Shift
    R8 = cinetica["wg_water_gas_shift"]["A"] * math.exp(-cinetica["wg_water_gas_shift"]["E"] / (Rg * T_g)) * (p_parcial["CO"] * p_parcial["H2O"] - (p_parcial["CO2"] * p_parcial["H2"]) / (0.02565 * math.exp(3966 / T_g)))
    
    # Quebra de Alcatrão (Cracking Tar)
    R9 = cinetica["p2_cracking_tar"]["A"] * math.exp(-cinetica["p2_cracking_tar"]["E"] / (Rg * T_g)) * max(0, rho_g["Tar"])
    
    taxas_reacao = {
        "SECAGEM": m,
        "Pirólise_1": Rp1,
        "G1_BOUDOUARD": R1,
        "G2_GASIFICAÇÃO_H2O": R2,
        "G3_GASIFICAÇÃO_H2": R3,
        "G4_REFORMAÇÃO_CH4": R4,
        "C1_OXIDAÇÃO_CH4": R5,
        "C2_OXIDAÇÃO_CO": R6,
        "C3_OXIDAÇÃO_H2": R7,
        "WGS": R8,
        "P2_TAR_CRACKING": R9
    }
    
    # 6. Balanços de Massa (Sólidos e Gases)
    M_char = materiais["Solidos"]["Char"]["MM"] / 1000.0
    M_H2O = materiais["Gases"]["H2O"]["MM"] / 1000.0
    st_p1 = cinetica["pirolise_p1"]["estequiometria"]
    
    R_B = st_p1["B"] * Rp1
    R_char = st_p1["Char"] * Rp1 - R1*M_char - R2*M_char - R3*M_char
    R_M = -m
    
    # Loop genérico para montar a taxa de cada gás baseada nas estequiometrias do JSON
    R_gases = {}
    for comp in gases_keys:
        taxa_pirolise = (st_p1.get(comp, 0) * Rp1 / (materiais["Gases"][comp]["MM"]/1000.0))
        # Adiciona a umidade evaporada apenas no H2O
        if comp == "H2O": taxa_pirolise += (m / M_H2O)
        
        taxa_reacoes = (
            cinetica["g1_boudouard"]["estequiometria"].get(comp, 0) * R1 +
            cinetica["g2_gasificacao_H2O"]["estequiometria"].get(comp, 0) * R2 +
            cinetica["g3_gasificacao_H2"]["estequiometria"].get(comp, 0) * R3 +
            cinetica["g4_reformacao_CH4"]["estequiometria"].get(comp, 0) * R4 +
            cinetica["c1_oxidacao_CH4"]["estequiometria"].get(comp, 0) * R5 +
            cinetica["c2_oxidacao_CO"]["estequiometria"].get(comp, 0) * R6 +
            cinetica["c3_oxidacao_H2"]["estequiometria"].get(comp, 0) * R7 +
            cinetica["wg_water_gas_shift"]["estequiometria"].get(comp, 0) * R8 +
            cinetica["p2_cracking_tar"]["estequiometria"].get(comp, 0) * (R9 / (materiais["Gases"][comp]["MM"]/1000.0))
        )
        R_gases[comp] = taxa_pirolise + taxa_reacoes

    # 7. Equações Diferenciais Ordinárias (EDOs)
    dydz = [0.0] * 16
    
    # Derivadas de velocidade (Sólido e Gás) e Pressão (Ergun)
    dUs_dz = -(1.0 / reator_params["rhochar0"]) * M_char * (R1 + R2 + R3)
    dydz[13] = dUs_dz
    
    dp_dz = -((150.0 * mu * ((1 - reator_params["eps"])**2) * U_g / ((reator_params["d_p"]**2) * (reator_params["eps"]**3))) +
              (1.75 * rhom_g * (1 - reator_params["eps"]) * (U_g**2) / (reator_params["d_p"] * (reator_params["eps"]**3))))
    dydz[12] = dp_dz

    # Derivadas de Massa (Sólidos)
    dydz[0] = (1.0 / U_s) * (R_B - rho_B * dUs_dz) if U_s > 1e-8 else 0.0
    dydz[1] = (1.0 / U_s) * (R_char - rho_char * dUs_dz) if U_s > 1e-8 else 0.0
    dydz[2] = (1.0 / U_s) * (R_M - rho_M * dUs_dz) if U_s > 1e-8 else 0.0
    drhoms_dz = dydz[0] + dydz[1] + dydz[2]
    
    # Derivadas de Massa (Gases)
    drhomg_dz = 0.0
    for i, comp in enumerate(gases_keys):
        dCg_dz = (1.0 / U_g) * (R_gases[comp] - C_g[i] * 0.0) # Simplificação antes de acoplar dUg_dz
        dydz[3 + i] = dCg_dz
        drhomg_dz += (materiais["Gases"][comp]["MM"]/1000.0) * dCg_dz
        
    dUg_dz = (-drhoms_dz * U_s - rhom_s * dUs_dz - drhomg_dz * U_g) / rhom_g
    dydz[11] = dUg_dz
    
    # Re-calcula derivadas dos gases com dUg_dz correto
    for i, comp in enumerate(gases_keys):
        dydz[3 + i] = (1.0 / U_g) * (R_gases[comp] - C_g[i] * dUg_dz)

    # ==========================================
    # 8. BALANÇO DE ENERGIA (AQUECIMENTO E REAÇÃO)
    # ==========================================
    T0 = 298.15
    M_B = materiais["Solidos"]["Biomassa"]["MM"] / 1000.0
    M_char = materiais["Solidos"]["Char"]["MM"] / 1000.0
    M_W = materiais["Solidos"]["Umidade"]["MM"] / 1000.0

    somaCscp = (rho_B/M_B)*cp_s["Biomassa"] + (rho_char/M_char)*cp_s["Char"] + (rho_M/M_W)*cp_s["Umidade"]
    somaCgcp = sum(C_g[i] * cp_g[comp] for i, comp in enumerate(gases_keys))
    
    somaCsdcpdz = (rho_char/M_char)*dcps_dT["Char"]
    somaCgdcpdz = sum(C_g[i] * dcp_g_dT[comp] for i, comp in enumerate(gases_keys))
    
    somacpdCsdz = cp_s["Biomassa"]*(dydz[0]/M_B) + cp_s["Char"]*(dydz[1]/M_char) + cp_s["Umidade"]*(dydz[2]/M_W)
    somacpdCgdz = sum(cp_g[comp]*dydz[3+i] for i, comp in enumerate(gases_keys))

    # Somatório do calor liberado/consumido pelas reações químicas
    somaRDH_s = Rp1*cinetica["pirolise_p1"]["dH"] + m*cinetica["secagem_m"]["dH"] + R1*cinetica["g1_boudouard"]["dH"] + R2*cinetica["g2_gasificacao_H2O"]["dH"] + R3*cinetica["g3_gasificacao_H2"]["dH"]
    somaRDH_g = R9*cinetica["p2_cracking_tar"]["dH"] + R5*cinetica["c1_oxidacao_CH4"]["dH"] + R6*cinetica["c2_oxidacao_CO"]["dH"] + R7*cinetica["c3_oxidacao_H2"]["dH"] + R8*cinetica["wg_water_gas_shift"]["dH"] + R4*cinetica["g4_reformacao_CH4"]["dH"]

    # dT_s / dz (Variação da Temperatura do Sólido)
    if rhom_s > 1e-4:
        denom_s = U_s * somaCscp + U_s * (T_s - T0) * somaCsdcpdz
        dydz[14] = (1.0 / denom_s) * (-(U_s*(T_s - T0)*somacpdCsdz + (T_s - T0)*somaCscp*dydz[13]) - somaRDH_s - transporte["Qsg"] - transporte["Qsw"])
    else:
        dydz[14] = 0.0

    # dT_g / dz (Variação da Temperatura do Gás)
    denom_g = U_g * somaCgcp + U_g * (T_g - T0) * somaCgdcpdz
    dydz[15] = (1.0 / denom_g) * (-(U_g*(T_g - T0)*somacpdCgdz + (T_g - T0)*somaCgcp*dydz[11]) - somaRDH_g + transporte["Qsg"] - transporte["Qgw"])

    return dydz, taxas_reacao