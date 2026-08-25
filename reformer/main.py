import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import reformer_odes

# ==========================================
# 1. LEITURA DOS DADOS DO GASEIFICADOR
# ==========================================
caminho_entrada = os.path.join(os.path.dirname(__file__), '../gasifier/', 'saida_gaseificador.json')
with open(caminho_entrada, 'r', encoding='utf-8') as f:
    dados_in = json.load(f)

# ==========================================
# 2. CARREGAMENTO DOS MATERIAIS E CONSTANTES
# ==========================================
# Reaproveitando o dicionário de materiais do gaseificador (se necessário pela cinética)
caminho_materiais = os.path.join(os.path.dirname(__file__), '..', 'gasifier', 'materials.json')
try:
    with open(caminho_materiais, 'r', encoding='utf-8') as f:
        materiais = json.load(f)
except FileNotFoundError:
    materiais = {} 

# Parâmetros Iniciais do Bi-Reformer (Estimativa para o modelo numérico)
reator_params = {
    "W_cat": 100.0,   # Massa de catalisador no leito (kg) - Variável de integração
    "D": 1.0,         # Diâmetro do reator (m)
    "porosidade": 0.4
}

# ==========================================
# 3. CONDIÇÕES INICIAIS (W = 0)
# ==========================================
saida_gas = dados_in["Saida_Gaseificador"]
vazoes = saida_gas["Vazoes_Molares_mol_s"]

T0 = saida_gas["Temperatura_K"]
P0 = saida_gas["Pressao_Pa"]

# Vazões Iniciais (mol/s) importadas do balanço exato do gaseificador
F_CO0 = vazoes.get("CO", 0.0)
F_H20 = vazoes.get("H2", 0.0)
F_H2O0 = vazoes.get("H2O", 0.0)
F_CO20 = vazoes.get("CO2", 0.0)
F_N20 = vazoes.get("N2", 0.0)

# --- VARIÁVEL DE MAKE-UP DE METANO ---
makeup_CH4 = 100.0  # Zerado para teste do cenário puramente acoplado
F_CH40 = vazoes.get("CH4", 0.0) + makeup_CH4

# Vetor de estado inicial: [F_CO, F_H2, F_CH4, F_H2O, F_CO2, F_N2, T, P]
y0 = [F_CO0, F_H20, F_CH40, F_H2O0, F_CO20, F_N20, T0, P0]

# ==========================================
# 4. INTEGRAÇÃO NUMÉRICA (SOLVE_IVP)
# ==========================================
print("Iniciando a integração do Bi-Reformer...")

solucao = solve_ivp(
    fun=reformer_odes,
    t_span=(0, reator_params["W_cat"]),
    y0=y0,
    args=(materiais, reator_params),
    method='RK45',    # Mantido em RK45 para avaliação preliminar do solver
    dense_output=True
)

print("Integração do Reformer concluída!")

W_vals = solucao.t
F_CO_vals = solucao.y[0]
F_H2_vals = solucao.y[1]
F_CH4_vals = solucao.y[2]
F_H2O_vals = solucao.y[3]
F_CO2_vals = solucao.y[4]
F_N2_vals = solucao.y[5]
T_vals = solucao.y[6]
P_vals = solucao.y[7]

# ==========================================
# 5. EXTRAÇÃO DE RESULTADOS E PLOTAGEM
# ==========================================
# Cálculo contínuo do Módulo Estequiométrico (S) ao longo do leito catalítico
S_vals = []
for h2, co2, co in zip(F_H2_vals, F_CO2_vals, F_CO_vals):
    denom = co + co2
    S = (h2 - co2) / denom if denom > 0 else 0
    S_vals.append(S)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Gráfico 1: Perfil Molar de Consumo e Geração
ax1.plot(W_vals, F_H2_vals, 'cyan', label='H₂', linewidth=2)
ax1.plot(W_vals, F_CO_vals, 'red', label='CO', linewidth=2)
ax1.plot(W_vals, F_CH4_vals, 'orange', label='CH₄', linewidth=2)
ax1.plot(W_vals, F_CO2_vals, 'magenta', label='CO₂', linewidth=2)
ax1.set_xlabel('Massa de Catalisador, W (kg)', fontsize=12)
ax1.set_ylabel('Vazão Molar (mol/s)', fontsize=12)
ax1.set_title('Perfil de Conversão no Bi-Reformer', fontsize=14)
ax1.legend()
ax1.grid(True)

# Gráfico 2: Adequação para o Reator de Metanol
ax2.plot(W_vals, S_vals, 'purple', linewidth=2)
ax2.axhline(y=2.0, color='green', linestyle='--', label='Alvo (S = 2.0)')
ax2.set_xlabel('Massa de Catalisador, W (kg)', fontsize=12)
ax2.set_ylabel('Módulo Estequiométrico (S)', fontsize=12)
ax2.set_title('Correção da Proporção do Syngas', fontsize=14)
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

# ==========================================
# 6. GERAÇÃO DO RELATÓRIO DE SAÍDA (JSON)
# ==========================================
relatorio_reformer = {
    "Entrada_Reformer": {
        "Temperatura_K": float(T0),
        "Pressao_Pa": float(P0),
        "Vazao_Makeup_CH4_mol_s": float(makeup_CH4)
    },
    "Saida_Reformer": {
        "Temperatura_K": float(T_vals[-1]),
        "Pressao_Pa": float(P_vals[-1]),
        "Modulo_S_Final": float(S_vals[-1]),
        "Vazoes_Molares_mol_s": {
            "CO": float(F_CO_vals[-1]),
            "H2": float(F_H2_vals[-1]),
            "CH4": float(F_CH4_vals[-1]),
            "H2O": float(F_H2O_vals[-1]),
            "CO2": float(F_CO2_vals[-1]),
            "N2": float(F_N2_vals[-1])
        }
    }
}

caminho_saida = os.path.join(os.path.dirname(__file__), 'saida_reformer.json')
with open(caminho_saida, 'w', encoding='utf-8') as f:
    json.dump(relatorio_reformer, f, indent=4)

print(f"\nRelatório do Reformer gerado com sucesso em: {caminho_saida}")