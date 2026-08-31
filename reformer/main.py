import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import cinetica_reformer

# ==========================================
# 1. LEITURA DOS DADOS DE SAÍDA DO GASEIFICADOR
# ==========================================

caminho_saida_gas = os.path.join(os.path.dirname(__file__), '..', 'gasifier', 'saida_gaseificador.json')

with open(caminho_saida_gas, 'r', encoding='utf-8') as f:
    dados_in = json.load(f)

saida = dados_in["Saida_Gaseificador"]
vazoes_in = dict(saida["Vazoes_Molares_mol_s"])  

# Removendo O2 e Tar (são zerados na saída do gaseificador)
vazoes_in.pop("O2", None)
vazoes_in.pop("Tar", None)

T_in_gasifier = saida["Temperatura_K"]
P_in = saida["Pressao_Pa"]

# ==========================================
# 2. CARREGAMENTO DOS MATERIAIS E CONSTANTES
# ==========================================
with open('./reformer/materiais.json', 'r') as materiais:
    componentes = json.load(materiais)


# ==========================================
# 3. PREPARAÇÃO DA ALIMENTAÇÃO DO REATOR WGS
# ==========================================
# Injeção de vapor d'água para corrigir a razão H2:CO do waste-derived syngas
#TODO: VERIFICAR A MELHOR FORMA DE FAZER A ADIÇÃO, QUERO FAZER UMA ANÁLISE DE SENSIBILIDADE APÓS ISSO
vazao_vapor_adicional = 5  # mol/s
vazoes_in["H2O"] = vazoes_in["H2O"] + vazao_vapor_adicional

# Injeção de CH4 adicional
vazao_CH4_adicional = 5 # mol/s
vazoes_in["CH4"] = vazoes_in["CH4"] + vazao_CH4_adicional

T_in = 503  # Típico para início de leito HTS (~350°C a 400°C)

Massa_Catalisador_Total = (0, 21) # kg de catalisador

# Vetor de estado [CO, H2, CH4, H2O, CO2, N2, T]
y0 = [
    vazoes_in.get("CO", 0), 
    vazoes_in.get("H2", 0), 
    vazoes_in.get("CH4", 0), 
    vazoes_in.get("H2O", 0), 
    vazoes_in.get("CO2", 0), 
    vazoes_in.get("N2", 0), 
    T_in,
    P_in,
    0
]

# ==========================================
# 4. DEFINIÇÃO DO SISTEMA DE EDOs
# ==========================================
def sistema_edos(w, y, componentes):
    # y representa as variáveis dependentes (vazão molar, por exemplo)
    # y é um vetor com as quantidades de cada componente no sistema
    condicoes_iniciais = {
        "CO": y[0],
        "H2": y[1],
        "CH4": y[2],
        "H2O": y[3],
        "CO2": y[4],
        "N2": y[5],
        "T": y[6],
        "P": y[7],
        "Q": y[8]
    }
    
    # Chama a função de cinética que calcula as taxas de geração
    taxa_de_geracao = cinetica_reformer(condicoes_iniciais, componentes)

    # Converte o dicionário de taxas de geração em uma lista de derivadas
    dy_dt = [
        taxa_de_geracao["CO"],
        taxa_de_geracao["H2"],
        taxa_de_geracao["CH4"],
        taxa_de_geracao["H2O"],
        taxa_de_geracao["CO2"],
        taxa_de_geracao["N2"],
        taxa_de_geracao["T"],
        taxa_de_geracao["P"],
        taxa_de_geracao["Q"]
    ]
    
    return dy_dt

# ==========================================
# 5. INTEGRAÇÃO NUMÉRICA
# ==========================================
print("Integrando reator WGS (correção H2:CO)...")
solucao = solve_ivp(
    sistema_edos,
    Massa_Catalisador_Total,
    y0,
    args=(componentes,),
    method='RK45'
)

massa_w = solucao.t
y_out = solucao.y

F_CO_out, F_H2_out, F_CH4_out, F_H2O_out, F_CO2_out, F_N2_out, T_out, P_out, Q_out = y_out[:, -1]

H2_CO_in = y0[1] / y0[0] if y0[0] > 0 else 0
H2_CO_out = F_H2_out / F_CO_out if F_CO_out > 0 else 0

# ==========================================
# 6. EXPORTAÇÃO DOS RESULTADOS (JSON)
# ==========================================
relatorio = {
    "Alimentacao": {
        "Vapor_Adicionado_mol_s": float(vazao_vapor_adicional),
        "Metano_Adicionado_mol_s": float(vazao_CH4_adicional),
        "Razao_H2_CO_Entrada": float(H2_CO_in),
        "Temperatura_Entrada_K": float(T_in)  # <-- CORREÇÃO 3: Removido o "_wgs"
    },
    "Saida": {
        "Vazoes_Molares": {
            "CO": float(F_CO_out), 
            "H2": float(F_H2_out), 
            "CH4": float(F_CH4_out),
            "H2O": float(F_H2O_out), 
            "CO2": float(F_CO2_out), 
            "N2": float(F_N2_out)
        },
        "Razao_H2_CO_Saida": float(H2_CO_out),
        "Temperatura_Saida_K": float(T_out),
        "Pressao_Saida_Pa": float(P_out)
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

ax1.plot(massa_w, y_out[0], label='CO', color='red', linewidth=2)
ax1.plot(massa_w, y_out[1], label='H₂', color='blue', linewidth=2)
ax1.plot(massa_w, y_out[2], label='CH₄', color='green', linewidth=2)
ax1.plot(massa_w, y_out[3], label='H₂O', color='cyan', linewidth=2)
ax1.plot(massa_w, y_out[4], label='CO₂', color='gray', linewidth=2)
ax1.set_xlabel('Massa de Catalisador (kg)')
ax1.set_ylabel('Vazão Molar (mol/s)')
ax1.legend(loc='center left')

ax2 = ax1.twinx()
ax2.plot(massa_w, y_out[6] - 273.15, label='Temperatura (°C)', color='orange', linestyle='--', linewidth=2)
ax2.set_ylabel('Temperatura (°C)')
ax2.legend(loc='center right')

plt.title('Perfil do Reformer (Ajuste de Gás de Síntese)')
plt.grid(True)
plt.show()