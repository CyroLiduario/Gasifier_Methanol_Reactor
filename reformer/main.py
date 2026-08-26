import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import cinetica_wgs

# ==========================================
# 1. LEITURA DOS DADOS DE SAÍDA DO GASEIFICADOR
# ==========================================
caminho_entrada = os.path.join(os.path.dirname(__file__), '../gasifier/', 'saida_gaseificador.json')

try:
    with open(caminho_entrada, 'r', encoding='utf-8') as f:
        dados_in = json.load(f)
    # Ajuste as chaves abaixo conforme a estrutura exata do seu 'saida_gaseificador.json'
    vazoes_in = dados_in.get("Vazoes_Molares_Saida", {"CO": 29.68, "H2": 19.07, "CO2": 7.87, "CH4": 0.8, "H2O": 3.72, "N2": 38.86})
    T_in_gasifier = dados_in.get("Temperatura_K", 1073.15)
    P_bar = dados_in.get("Pressao_bar", 25.0)
except FileNotFoundError:
    print("Aviso: 'saida_gaseificador.json' não encontrado. Utilizando dados padrão do syngas.")
    vazoes_in = {"CO": 29.68, "H2": 19.07, "CO2": 7.87, "CH4": 0.8, "H2O": 3.72, "N2": 38.86}
    T_in_gasifier = 950.0 + 273.15
    P_bar = 25.0

# ==========================================
# 2. CARREGAMENTO DOS MATERIAIS E CONSTANTES
# ==========================================
caminho_materiais = os.path.join(os.path.dirname(__file__), '..', 'gasifier', 'materials.json')
try:
    with open(caminho_materiais, 'r', encoding='utf-8') as f:
        materiais = json.load(f)
except FileNotFoundError:
    materiais = {} 

# ==========================================
# 3. PREPARAÇÃO DA ALIMENTAÇÃO DO REATOR WGS
# ==========================================
# Injeção de vapor d'água para corrigir a razão H2:CO do waste-derived syngas
vazao_vapor_adicional = 30.0  # mol/s (Ajuste este valor no controle do seu processo)
vazoes_in["H2O"] += vazao_vapor_adicional

# A temperatura de entrada no reator HTS (High-Temperature Shift) costuma ser menor que a do gaseificador
T_in_wgs = 350.0 + 273.15  # Típico para início de leito HTS (~350°C a 400°C)

Massa_Catalisador_Total = 8000.0 # kg de catalisador

# Vetor de estado [CO, H2, CH4, H2O, CO2, N2, T]
y0 = [
    vazoes_in.get("CO", 0), 
    vazoes_in.get("H2", 0), 
    vazoes_in.get("CH4", 0), 
    vazoes_in.get("H2O", 0), 
    vazoes_in.get("CO2", 0), 
    vazoes_in.get("N2", 0), 
    T_in_wgs
]

# ==========================================
# 4. DEFINIÇÃO DO SISTEMA DE EDOs
# ==========================================
def balanco_reator_wgs(w, y):
    F_CO, F_H2, F_CH4, F_H2O, F_CO2, F_N2, T = y
    
    # Obter taxas da cinética
    taxa = cinetica_wgs(y, P_bar)
    
    # Capacidades caloríficas simplificadas (J/mol.K)
    Cp = {"CO": 30.0, "H2": 29.0, "CH4": 50.0, "H2O": 36.0, "CO2": 45.0, "N2": 29.0}
    sum_FCp = (F_CO*Cp["CO"] + F_H2*Cp["H2"] + F_CH4*Cp["CH4"] + 
               F_H2O*Cp["H2O"] + F_CO2*Cp["CO2"] + F_N2*Cp["N2"])
    
    # Balanço de energia adiabático (Q = 0)
    dT_dw = taxa["calor"] / sum_FCp
    
    return [
        taxa["CO"],
        taxa["H2"],
        taxa["CH4"],
        taxa["H2O"],
        taxa["CO2"],
        taxa["N2"],
        dT_dw
    ]

# ==========================================
# 5. INTEGRAÇÃO NUMÉRICA
# ==========================================
print("Integrando reator WGS (correção H2:CO)...")
solucao = solve_ivp(
    balanco_reator_wgs,
    (0, Massa_Catalisador_Total),
    y0,
    method='BDF',
    dense_output=True
)

massa_w = solucao.t
y_out = solucao.y

F_CO_out, F_H2_out, F_CH4_out, F_H2O_out, F_CO2_out, F_N2_out, T_out = y_out[:, -1]
H2_CO_in = y0[1] / y0[0]
H2_CO_out = F_H2_out / F_CO_out

# ==========================================
# 6. EXPORTAÇÃO DOS RESULTADOS (JSON)
# ==========================================
relatorio = {
    "Alimentacao": {
        "Vapor_Adicionado_mol_s": vazao_vapor_adicional,
        "Razao_H2_CO_Entrada": H2_CO_in,
        "Temperatura_Entrada_K": T_in_wgs
    },
    "Saida": {
        "Vazoes_Molares": {
            "CO": F_CO_out, "H2": F_H2_out, "CH4": F_CH4_out,
            "H2O": F_H2O_out, "CO2": F_CO2_out, "N2": F_N2_out
        },
        "Razao_H2_CO_Saida": H2_CO_out,
        "Temperatura_Saida_K": T_out
    }
}

caminho_saida = os.path.join(os.path.dirname(__file__), 'saida_wgs.json')
with open(caminho_saida, 'w', encoding='utf-8') as f:
    json.dump(relatorio, f, indent=4)

print(f"Sucesso! Relatório salvo em: {caminho_saida}")
print(f"Razão H2:CO ajustada de {H2_CO_in:.2f} para {H2_CO_out:.2f}")

# ==========================================
# 7. VISUALIZAÇÃO DOS PERFIS
# ==========================================
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(massa_w, y_out[0], label='CO', color='red')
ax1.plot(massa_w, y_out[1], label='H₂', color='blue')
ax1.plot(massa_w, y_out[3], label='H₂O', color='cyan')
ax1.plot(massa_w, y_out[4], label='CO₂', color='gray')
ax1.set_xlabel('Massa de Catalisador (kg)')
ax1.set_ylabel('Vazão Molar (mol/s)')
ax1.legend(loc='center left')

ax2 = ax1.twinx()
ax2.plot(massa_w, y_out[6] - 273.15, label='Temperatura (°C)', color='orange', linestyle='--')
ax2.set_ylabel('Temperatura (°C)')
ax2.legend(loc='center right')

plt.title('Perfil do Reator WGS (Alimentado pelo Gaseificador)')
plt.grid(True)
plt.show()