# Realiza o cálculo de um reator de bi-reforma catalítica 
# conforme proposto por Ramos et al. (2024) e García-Vargas et al. (2015)
'''
    Catalyst: Ni-Mg/β-SiC
    Dry Reforming: CH4 + CO2 ↔ 2 CO + 2 H2
    Steam Reforming: CH4 + H2O ↔ CO+ 3 H2
    Water-gas shift: CO + H2O ↔ CO2 + H2
    Reverse Methanation: CH4 + 2 H2O ↔ CO2 + 4 H2 - desconsiderada no estudo de Ramos et al. (2024)
'''

#TODO: CONFIRMAR DADOS DE DENSIDADE DO CATALIZADOR

import os
import json
from math import pi
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

# Removendo Tar (são zerados na saída do gaseificador)
vazoes_in.pop("Tar", None)

T_in = saida["Temperatura_K"]
P_in = saida["Pressao_Pa"]

# ==========================================
# 2. CARREGAMENTO DOS MATERIAIS E CONSTANTES
# ==========================================
with open('./reformer_acquarola/materiais.json', 'r') as materiais:
    componentes = json.load(materiais)

# ==========================================
# 3. PREPARAÇÃO DA ALIMENTAÇÃO DO REATOR DE REFORMA
# ==========================================
# Injeção de vapor d'água para corrigir a razão H2:CO do waste-derived syngas
#TODO: VERIFICAR A MELHOR FORMA DE FAZER A ADIÇÃO, QUERO FAZER UMA ANÁLISE DE SENSIBILIDADE APÓS ISSO
vazao_vapor_adicional = 216.44  # mol/s
vazoes_in["H2O"] = vazoes_in["H2O"] + vazao_vapor_adicional

# Injeção de CH4 adicional
vazao_CH4_adicional = 216.44 # mol/s
vazoes_in["CH4"] = vazoes_in["CH4"] + vazao_CH4_adicional

multitubular = False
if multitubular:
    numero_de_tubos = 3000
    diametro_reator = 0.0762
    vazoes_in = {componente: vazao_molar_comp/numero_de_tubos for componente, vazao_molar_comp in vazoes_in.items()}
else:
    diametro_reator = 3
    numero_de_tubos = 1

Massa_Catalisador_Total = (0, 1000) # kg de catalisador


# Vetor de estado [CO, H2, CH4, H2O, CO2, N2, O2, T]
y0 = [
    vazoes_in.get("CO", 0), 
    vazoes_in.get("H2", 0), 
    vazoes_in.get("CH4", 0), 
    vazoes_in.get("H2O", 0), 
    vazoes_in.get("CO2", 0), 
    vazoes_in.get("N2", 0), 
    vazoes_in.get("O2", 0),
    T_in,
    P_in,
    0
]

# ==========================================
# 4. DEFINIÇÃO DO SISTEMA DE EDOs
# ==========================================
def sistema_edos(w, y, componentes, diametro_reator):
    # y representa as variáveis dependentes (vazão molar, por exemplo)
    # y é um vetor com as quantidades de cada componente no sistema
    condicoes_iniciais = {
        "CO": y[0],
        "H2": y[1],
        "CH4": y[2],
        "H2O": y[3],
        "CO2": y[4],
        "N2": y[5],
        "O2": y[6],
        "T": y[7],
        "P": y[8],
        "Q": y[9]
    }
    
    # Chama a função de cinética que calcula as taxas de geração
    taxa_de_geracao, _ = cinetica_reformer(condicoes_iniciais, componentes, diametro_reator)

    # Converte o dicionário de taxas de geração em uma lista de derivadas
    dy_dt = [
        taxa_de_geracao["CO"],
        taxa_de_geracao["H2"],
        taxa_de_geracao["CH4"],
        taxa_de_geracao["H2O"],
        taxa_de_geracao["CO2"],
        taxa_de_geracao["N2"],
        taxa_de_geracao["O2"],
        taxa_de_geracao["T"],
        taxa_de_geracao["P"],
        taxa_de_geracao["Q"]
    ]
    
    return dy_dt

# ==========================================
# 5. INTEGRAÇÃO NUMÉRICA
# ==========================================
print("Integrando reator de reforma (correção H2:CO)...")
solucao = solve_ivp(
    sistema_edos,
    Massa_Catalisador_Total,
    y0,
    args=(componentes,diametro_reator),
    method='RK45'
)

massa_w = solucao.t
y_out = solucao.y

F_CO_out, F_H2_out, F_CH4_out, F_H2O_out, F_CO2_out, F_N2_out, F_O2_out, T_out, P_out, Q_out = y_out[:, -1]*numero_de_tubos

H2_CO_in = y0[1] / y0[0] if y0[0] > 0 else 0
H2_CO_out = F_H2_out / F_CO_out if F_CO_out > 0 else 0


# ==========================================
# 6. PERFIS DAS TAXAS DE REAÇÃO
# ==========================================

taxas_drm = []
taxas_srm = []
taxas_wgs = []

for i in range(len(massa_w)):

    condicoes_iniciais = {
        "CO": y_out[0, i],
        "H2": y_out[1, i],
        "CH4": y_out[2, i],
        "H2O": y_out[3, i],
        "CO2": y_out[4, i],
        "N2": y_out[5, i],
        "O2": y_out[6, i],
        "T": y_out[7, i],
        "P": y_out[8, i],
        "Q": y_out[9, i]
    }

    _, taxas_reacao = cinetica_reformer(
        condicoes_iniciais,
        componentes,
        diametro_reator
    )

    taxas_drm.append(taxas_reacao["DRM"])
    taxas_srm.append(taxas_reacao["SRM"])
    taxas_wgs.append(taxas_reacao["WGS"])

# ==========================================
# 7. EXPORTAÇÃO DOS RESULTADOS (JSON)
# ==========================================
relatorio = {
    "Alimentacao": {
        "Vapor_Adicionado_mol_s": float(vazao_vapor_adicional),
        "Metano_Adicionado_mol_s": float(vazao_CH4_adicional),
        "Razao_H2_CO_Entrada": float(H2_CO_in),
        "Temperatura_Entrada_K": float(T_in),
        "Vazoes_Molares": vazoes_in
    },
    "Saida": {
        "Vazoes_Molares": {
            "CO": float(F_CO_out), 
            "H2": float(F_H2_out), 
            "CH4": float(F_CH4_out),
            "H2O": float(F_H2O_out), 
            "CO2": float(F_CO2_out), 
            "N2": float(F_N2_out),
            "O2": float(F_O2_out)
        },
        "Razao_H2_CO_Saida": float(H2_CO_out),
        "Temperatura_Saida_K": float(T_out),
        "Pressao_Saida_Pa": float(P_out)
    }
}

caminho_saida = os.path.join(os.path.dirname(__file__), 'saida_reformer.json')
with open(caminho_saida, 'w', encoding='utf-8') as f:
    json.dump(relatorio, f, indent=4)

print(f"Sucesso! Relatório salvo em: {caminho_saida}")
print(f"Razão H2:CO ajustada de {H2_CO_in:.2f} para {H2_CO_out:.2f}")

# ==========================================
# 8. VISUALIZAÇÃO DOS PERFIS
# ==========================================

# ===============================
#   GRÁFICOS DE TAXA DE REAÇÃO
# ===============================
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# ==========================
# DRY REFORMING
# ==========================
axes[0].plot(
    massa_w,
    taxas_drm,
    linewidth=1.5
)

axes[0].set_xlabel("Massa de catalisador [kg]")
axes[0].set_ylabel("Taxa [mol/(kg·s)]")
axes[0].set_title("Dry Reforming")
axes[0].grid(True, which="both")


# ==========================
# STEAM REFORMING
# ==========================
axes[1].plot(
    massa_w,
    taxas_srm,
    linewidth=1.5
)

axes[1].set_xlabel("Massa de catalisador [kg]")
axes[1].set_ylabel("Taxa [mol/(kg·s)]")
axes[1].set_title("Steam Reforming")
axes[1].grid(True, which="both")


# ==========================
# WATER-GAS SHIFT
# ==========================
axes[2].plot(
    massa_w,
    taxas_wgs,
    linewidth=1.5
)

axes[2].set_xlabel("Massa de catalisador [kg]")
axes[2].set_ylabel("Taxa [mol/(kg·s)]")
axes[2].set_title("Water-Gas Shift")
axes[2].grid(True, which="both")


fig.suptitle("Perfil das taxas de reação ao longo do leito")
plt.tight_layout()


# ==========================
#   GRÁFICO DE COMPOSIÇÃO
# ==========================
fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(massa_w, y_out[0], label='CO', color='red', linewidth=2)
ax1.plot(massa_w, y_out[1], label='H₂', color='blue', linewidth=2)
ax1.plot(massa_w, y_out[2], label='CH₄', color='green', linewidth=2)
ax1.plot(massa_w, y_out[3], label='H₂O', color='cyan', linewidth=2)
ax1.plot(massa_w, y_out[4], label='CO₂', color='gray', linewidth=2)
ax1.set_xlabel('Massa de Catalisador (kg)')
ax1.set_ylabel('Vazão Molar (mol/s)')
ax1.legend(loc='center left')

plt.title('Perfil do Reformer (Ajuste de Gás de Síntese)')
plt.grid(True)

# ==========================
#   PERFIL DE PRESSÃO
# ==========================
pressao = solucao.y[7]
# Plotando os resultados
plt.figure(figsize=(10, 6))
plt.plot(massa_w, pressao, label="P [Pa]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Pressão [Pa]', fontsize=12)
plt.title('Perfil de Pressão', fontsize=14)
plt.legend()
plt.grid(True)

# ==========================
#   PERFIL DE TEMPERATURA
# ==========================

temperatura = solucao.y[6]
# Plotando os resultados
plt.figure(figsize=(10, 6))

plt.plot(massa_w, temperatura, label="T [K]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Temperatura [K]', fontsize=12)
plt.title('Perfil de Temperatura', fontsize=14)
plt.legend()
plt.grid(True)

# ==========================
#   PERFIL DE CALOR
# ==========================

calor = solucao.y[8]
# Plotando os resultados
plt.figure(figsize=(10, 6))

plt.plot(massa_w, calor, label="Q [J/s]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Calor [J]', fontsize=12)
plt.title('Perfil de Calor', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

# ==========================
#   TAMANHO DO REATOR
# ==========================

# ALTURA
area_secao_transversal = pi*(diametro_reator/2)**2  # m² - tubo de 2 in

porosidade = 0.6 # suposta - Salles (2013)
densidade_particula_catalisador = 3950  # kg_cat/m³ - Salles (2013)
diametro_particula = 0.001  # m - GARCÍA-VARGAS et al. (2015)

densidade_do_leito = densidade_particula_catalisador * (1-porosidade)
volume = max(massa_w)/densidade_do_leito
altura_reator = volume/area_secao_transversal

print(f"Altura do Reator = {altura_reator:.2f} m")

# DIAMETRO

if multitubular:
    # Usando um passe triangular
    diametro_externo = 0.06 # m - confirmar, considerei tubos de 2 in e e schedule 40 
    pitch = 1.25*diametro_externo # distância entre os centros dos tubos considerando arranjo triangular 
    diametro_casco = pitch*((numero_de_tubos/0.315)**0.5)

else:
    diametro_casco = diametro_reator

print(f"Diâmetro do Reator = {diametro_casco:.2f} m")