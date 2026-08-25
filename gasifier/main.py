import os
import json
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import gasifier_odes

# ==========================================
# 1. CARREGAMENTO DOS MATERIAIS E CONSTANTES
# ==========================================
with open('./gasifier/materials.json', 'r') as f:
    materiais = json.load(f)

# Constantes Fundamentais
Rg = 8.314 # J/(mol.K)
P0 = 3.0 * 101325.0 # Pressão inicial: 3 atm (Pa)
T0 = 300.0 # Temperatura inicial da biomassa e gás (K)

# Parâmetros Estruturais e Operacionais do Reator (Modelo Cyro Industrial)
reator_params = {
    "L": 12.0,            # Comprimento do reator (m)
    "D": 2.0,             # Diâmetro do reator (m)
    "eps": 0.5,           # Porosidade do leito
    "phi": 2.0 / 3.0,     # Fator de forma
    "em": 0.9,            # Emissividade da biomassa
    "d_p": 0.042544,      # Diâmetro equivalente da partícula (m)
    "rhochar0": 152.5,    # Densidade inicial do carvão vegetal (kg/m³)
    "sigma": 5.67051e-8,  # Constante de Boltzmann (W/m².K⁴)
    "T_w": 1050.0         # Temperatura da parede do reator (K)
}

# ==========================================
# 2. CONDIÇÕES INICIAIS (z = 0)
# ==========================================
area_secao = (math.pi * reator_params["D"]**2) / 4.0

# --- Fase Sólida ---
F_s = 101104.0 / 5.0 # Vazão mássica de biomassa (kg/h)
YM = 0.0872 # Umidade inicial (8.72%)
rhom_s0 = 119.0 # Densidade da biomassa (kg/m³)

U_s0 = (F_s / 3600.0) / (rhom_s0 * area_secao) # m/s
rho_B0 = rhom_s0 * (1.0 - YM)
rho_char0 = 0.0 # Sem carvão no topo
rho_M0 = rhom_s0 * YM

# --- Fase Gasosa (Ar de Combustão) ---
ER = 0.18 # Razão de equivalência

# Defina a fração molar de vapor d'água que você quer na mistura de entrada
#fracao_vapor = 0 # Exemplo: 15% de vapor e 85% de ar
#fracao_ar = 1.0 - fracao_vapor
#fracao_O2 = 0.21 * fracao_ar
#fracao_N2 = 0.79 * fracao_ar

# A ordem do vetor é: [H2O, O2, CO, CO2, H2, CH4, N2, Tar]
#X_g0 = [fracao_vapor, fracao_O2, 0.0, 0.0, 0.0, 0.0, fracao_N2, 0.0]
X_g0 = [0.0, 0.21, 0.0, 0.0, 0.0, 0.0, 0.79, 0.0]

rhom_g0_ar = 101325.0 / (287.058 * 300.0) # Densidade padrão do ar (kg/m³)

# Cálculo EXATO da estequiometria do Bagaço (Dados Cyro)
C, H, O, N, S = 43.02, 5.85, 44.34, 0.56, 0.89
w_bag = 1.0
x_bag = (H * 12.0) / (C * 1.0)
y_bag = (O * 12.0) / (C * 16.0)
z_bag = (N * 12.0) / (C * 14.0)
k_bag = (S * 12.0) / (C * 32.0)

mcomb_bag = (w_bag*12.0 + x_bag*1.0 + y_bag*16.0 + z_bag*14.0 + k_bag*32.0) * 1e-3 # kg
var_bag = (w_bag + x_bag/4.0 - y_bag/2.0) * (2*16.0 + 3.76*2*14.0) * 1e-3 / rhom_g0_ar # m³

ER_stoic = var_bag / mcomb_bag # ~4.59 m³ ar / kg biomassa

Q_g = ER * ER_stoic * F_s # Vazão volumétrica de ar (m³/h)
F_g = Q_g * rhom_g0_ar # Vazão mássica de ar (kg/h)

U_g0 = (F_g / 3600.0) / (rhom_g0_ar * area_secao) # Velocidade inicial do gás (m/s)

# Concentrações Molares do Gás Inicial (mol/m³)
Mm_g0 = (materiais["Gases"]["O2"]["MM"] * 0.21 + materiais["Gases"]["N2"]["MM"] * 0.79) / 1000.0
rhom_g0_real = (P0 * Mm_g0) / (Rg * T0)
Cm_g0 = rhom_g0_real / Mm_g0
C_g0 = [Cm_g0 * x for x in X_g0]

# Vetor de Estado Inicial: [rho_B, rho_char, rho_M, C_g_H2O...Tar, U_g, P, U_s, T_s, T_g]
y0 = [rho_B0, rho_char0, rho_M0] + C_g0 + [U_g0, P0, U_s0, T0, T0]

# ==========================================
# 3. INTEGRAÇÃO NUMÉRICA (SOLVE_IVP)
# ==========================================
print("Iniciando a integração do gaseificador industrial...")

solucao = solve_ivp(
    fun=gasifier_odes,
    t_span=(0, reator_params["L"]),
    y0=y0,
    args=(materiais, reator_params),
    method='RK45',    # <-- Método idêntico ao ode23 do MATLAB
    max_step=0.05,    # <-- Força o Python a ler o reator a cada 5 cm no máximo
    dense_output=True
)

z_vals = solucao.t
y_vals = solucao.y

print("Integração concluída com sucesso!")

# ==========================================
# 4. EXTRAÇÃO DE RESULTADOS E PLOTAGEM
# ==========================================
rho_solidos = y_vals[0:3]
C_gases = y_vals[3:11]
Ug_vals = y_vals[11]
P_vals = y_vals[12]
Us_vals = y_vals[13]
Ts_vals = y_vals[14]
Tg_vals = y_vals[15]

gases_labels = ["H₂O", "O₂", "CO", "CO₂", "H₂", "CH₄", "N₂", "Alcatrão"]
cores = ['blue', 'green', 'red', 'magenta', 'cyan', 'orange', 'gray', 'black']

# Gráfico 1: Concentração dos Gases
plt.figure(figsize=(10, 6))
for i in range(8):
    plt.plot(z_vals, C_gases[i], label=gases_labels[i], color=cores[i], linewidth=2)
plt.xlabel('Comprimento do Reator, z (m)', fontsize=12)
plt.ylabel('Concentração (mol/m³)', fontsize=12)
plt.title('Perfil de Concentração dos Gases', fontsize=14)
plt.legend()
plt.grid(True)

# Gráfico 2: Densidade dos Sólidos
plt.figure(figsize=(10, 6))
plt.plot(z_vals, rho_solidos[0], label="Biomassa", linewidth=2)
plt.plot(z_vals, rho_solidos[1], label="Carvão (Char)", linewidth=2)
plt.plot(z_vals, rho_solidos[2], label="Umidade", linewidth=2)
plt.xlabel('Comprimento do Reator, z (m)', fontsize=12)
plt.ylabel('Densidade (kg/m³)', fontsize=12)
plt.title('Perfil de Consumo da Fase Sólida', fontsize=14)
plt.legend()
plt.grid(True)

# Gráfico 3: Perfil Térmico
plt.figure(figsize=(10, 6))
plt.plot(z_vals, Tg_vals, label="Temperatura do Gás (Tg)", color='red', linewidth=2)
plt.plot(z_vals, Ts_vals, label="Temperatura do Sólido (Ts)", color='blue', linestyle='--', linewidth=2)
plt.xlabel('Comprimento do Reator, z (m)', fontsize=12)
plt.ylabel('Temperatura (K)', fontsize=12)
plt.title('Perfil de Temperaturas', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

# ==========================================
# 5. CÁLCULO MECÂNICO DA ESPESSURA (EXTRA CYRO)
# ==========================================
Pi = (P0 - 101325.0) / 101325.0 # Pressão manométrica em atm
Pd = math.exp(0.39303 + (0.92452 * math.log(Pi)) + (0.0015655 * (math.log(Pi)**2)))
sigma_aco = 910.0 # Tensão admissível em bar
ts = (1000.0 * Pd * reator_params["D"]) / ((2.0 * sigma_aco) - (1.2 * Pd))

print(f"Cálculo Estrutural:")
print(f"Espessura mínima requerida para a parede de aço = {ts:.2f} mm")

# ==========================================
# 6. GERAÇÃO DO RELATÓRIO DE SAÍDA (JSON)
# ==========================================

# 1. Cálculo das vazões molares de saída (mol/s)
# Q_out (m³/s) = U_g_out (m/s) * area_secao (m²)
U_g_out = Ug_vals[-1]
Q_out = U_g_out * area_secao 

chaves_gases = ["H2O", "O2", "CO", "CO2", "H2", "CH4", "N2", "Tar"]
vazoes_saida_mol_s = {}

for i, comp in enumerate(chaves_gases):
    # C_gases[i][-1] está em mol/m³. Multiplicando por m³/s temos mol/s
    vazoes_saida_mol_s[comp] = float(C_gases[i][-1] * Q_out)

# 2. Cálculo das vazões molares de entrada de gás (mol/s) para o relatório
# O índice 0 é H2O, 1 é O2, e 6 é N2 no seu array X_g0 / C_g0
Q_in = U_g0 * area_secao
vazao_vapor_in_mol_s = float(C_g0[0] * Q_in)
vazao_ar_in_mol_s = float((C_g0[1] + C_g0[6]) * Q_in) # Somando O2 e N2

# 3. Estruturação do Dicionário de Relatório
relatorio_integracao = {
    "Entrada_Gaseificador": {
        "Temperatura_K": float(T0),
        "Pressao_Pa": float(P0),
        "Vazao_Biomassa_kg_h": float(F_s),
        "Vazao_Ar_kg_h": float(F_g),
        "Razao_Equivalencia_ER": float(ER),
        "Vazao_Vapor_mol_s": vazao_vapor_in_mol_s,
        "Vazao_Ar_mol_s": vazao_ar_in_mol_s
    },
    "Saida_Gaseificador": {
        "Temperatura_K": float(Tg_vals[-1]),
        "Pressao_Pa": float(P_vals[-1]),
        "Vazoes_Molares_mol_s": vazoes_saida_mol_s
    }
}

# 4. Salvando o arquivo JSON na raiz do projeto
caminho_arquivo = os.path.join(os.path.dirname(__file__), 'saida_gaseificador.json')

with open(caminho_arquivo, 'w', encoding='utf-8') as f:
    json.dump(relatorio_integracao, f, indent=4)

print(f"\nRelatório de integração gerado com sucesso em: {caminho_arquivo}")