import math
from thermodynamic import capacidade_calorifica_gases, derivada_capacidade_calorifica_gases

def reformer_odes(W, y, materiais, reator_params):
    # y = [F_CO, F_H2, F_CH4, F_H2O, F_CO2, F_N2, T, P] (Vazões molares, Temp, Pressão)
    F_CO, F_H2, F_CH4, F_H2O, F_CO2, F_N2, T, P = y
    
    # Constantes
    R = 8.314 # J/mol.K
    
    # 1. Frações Molares e Pressões Parciais (kPa)
    F_total = F_CO + F_H2 + F_CH4 + F_H2O + F_CO2 + F_N2
    P_kPa = P / 1000.0 # O modelo de Ramos usa kPa ou bar dependendo do K_eq, verifique no artigo
    
    p_CO  = (F_CO / F_total) * P_kPa
    p_H2  = (F_H2 / F_total) * P_kPa
    p_CH4 = (F_CH4 / F_total) * P_kPa
    p_H2O = (F_H2O / F_total) * P_kPa
    p_CO2 = (F_CO2 / F_total) * P_kPa
    
    # 2. Constantes Cinéticas (k) e de Equilíbrio (K_eq) de Ramos et al.
    # Lembre-se de verificar as unidades da energia de ativação (se kJ/mol, multiplique R por 10^-3)
    k_DR = 70.99e-3 * math.exp(-77.82 / (0.008314 * T))
    K_DR = 6.780e18 * math.exp(-31230 / T)
    
    k_SR = 85.77e-9 * math.exp(-74.72 / (0.008314 * T))
    K_SR = 1.198e17 * math.exp(-26830 / T)
    
    k_WGS = 149.22e-3 * math.exp(-34.26 / (0.008314 * T))
    K_WGS = 10**((2073 / T) - 2.029)
    
    # 3. Taxas de Reação (r = k * Força Motriz)
    # Proteção contra divisão por zero e pressões negativas
    r_DR = k_DR * p_CH4 * (1 - ((p_CO**2 * p_H2**2) / max(1e-10, p_CH4 * p_CO2 * K_DR)))
    r_SR = k_SR * p_CH4 * (1 - ((p_CO * p_H2**3) / max(1e-10, p_CH4 * p_H2O * K_SR)))
    r_WGS = k_WGS * ( (p_CO * p_H2O / max(1e-10, p_H2)) - (p_CO2 / K_WGS) )
    
    # 4. Balanço de Massa (dF/dW)
    dF_CO_dW  = 2*r_DR + r_SR - r_WGS
    dF_H2_dW  = 2*r_DR + 3*r_SR + r_WGS
    dF_CH4_dW = -r_DR - r_SR
    dF_H2O_dW = -r_SR - r_WGS
    dF_CO2_dW = -r_DR + r_WGS
    dF_N2_dW  = 0.0
    
    # 5. Balanço de Energia (dT/dW) e Perda de Carga (dP/dW)
    # (Inserir os cálculos baseados na entalpia de reação e Cp, similar ao metanol)
    dT_dW = 0.0 # Substituir pelo balanço real
    dP_dW = 0.0 # Substituir pela equação de Ergun
    
    return [dF_CO_dW, dF_H2_dW, dF_CH4_dW, dF_H2O_dW, dF_CO2_dW, dF_N2_dW, dT_dW, dP_dW]