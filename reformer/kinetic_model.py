import numpy as np

def cinetica_wgs(y, P_total_bar):
    """
    Avalia as taxas de reação LHHW para o reator WGS.
    """
    F_CO, F_H2, F_CH4, F_H2O, F_CO2, F_N2, T = y
    
    F_total = F_CO + F_H2 + F_CH4 + F_H2O + F_CO2 + F_N2
    
    # Pressões parciais (bar)
    p = {
        "CO":  (F_CO / F_total) * P_total_bar,
        "H2":  max(1e-6, (F_H2 / F_total) * P_total_bar), # Evita divisão por zero
        "CH4": (F_CH4 / F_total) * P_total_bar,
        "H2O": (F_H2O / F_total) * P_total_bar,
        "CO2": (F_CO2 / F_total) * P_total_bar
    }
    
    Rg = 8.314 # J/(mol.K)
    
    # 1. Constante Cinética (k) e de Equilíbrio (K_eq) para WGS
    # CO + H2O <=> CO2 + H2
    k_wgs = 1.955e6 * np.exp(-67130 / (Rg * T))
    K_eq_wgs = np.exp(4400 / T - 4.036)
    
    # 2. Constantes de Adsorção (K)
    K_CO  = 8.23e-5 * np.exp(70650 / (Rg * T))
    K_H2  = 6.12e-9 * np.exp(82900 / (Rg * T))
    K_CH4 = 6.65e-4 * np.exp(38280 / (Rg * T))
    K_H2O = 1.77e5  * np.exp(-88680 / (Rg * T))
    
    # 3. Denominador da Força Motriz (Adsorção Competitiva em Sítio Ativo)
    DEN = 1.0 + K_CO*p["CO"] + K_H2*p["H2"] + K_CH4*p["CH4"] + (K_H2O * p["H2O"] / p["H2"])
    
    # 4. Taxa de Reação WGS (mol / kg_cat . s)
    r_wgs = (k_wgs / p["H2"]) * (p["CO"]*p["H2O"] - (p["H2"] * p["CO2"] / K_eq_wgs)) / (DEN**2)
    
    # 5. Balanço Estequiométrico das Espécies (mol / kg_cat . s)
    taxa_geracao = {
        "CO":  -r_wgs,
        "H2O": -r_wgs,
        "CO2":  r_wgs,
        "H2":   r_wgs,
        "CH4":  0.0,
        "N2":   0.0 
    }
    
    # Entalpia da reação WGS (~ -41.2 kJ/mol a -41.98 kJ/mol)
    dH_wgs = -41150.0  # J/mol (Reação exotérmica)
    
    taxa_geracao["calor"] = -(r_wgs * dH_wgs) # Calor liberado: J / (kg_cat . s)
    
    return taxa_geracao