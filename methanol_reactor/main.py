import os
import json
import matplotlib.pyplot as plt
from math import pi
from kinetic_model import cinetica
from scipy.integrate import solve_ivp

with open('./materiais.json', 'r') as materiais:
    componentes = json.load(materiais)

caminho_saida_gas = os.path.join(os.path.dirname(__file__), '..', 'reformer', 'saida_reformer.json')
#caminho_saida_gas = os.path.join(os.path.dirname(__file__), '..', 'reformer_optimization', 'saida_reformer.json')

with open(caminho_saida_gas, 'r', encoding='utf-8') as f:
    dados_gaseificador = json.load(f)

vazoes_brutas = dados_gaseificador["Saida"]["Vazoes_Molares"]

# 3. Mapeamento das correntes (mol/s globais)
vazao_molar = {
    "CO": vazoes_brutas.get("CO"),
    "H2": vazoes_brutas.get("H2"),
    "CH3OH": 0.0,  # O gaseificador não produz metanol
    "CH4": vazoes_brutas.get("CH4"),
    "H2O": vazoes_brutas.get("H2O"),
    "CO2": vazoes_brutas.get("CO2"),
    "N2": vazoes_brutas.get("N2")
}
wt_agua = vazao_molar["H2O"]/sum(vazao_molar.values())
print(wt_agua)
# vazao por tubo
numero_de_tubos = 300
vazao_molar = {componente: vazao_molar_comp/numero_de_tubos for componente, vazao_molar_comp in vazao_molar.items()}
#print(vazao_molar)

vazao_massica = {} # g/s
for componente in componentes:
    vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

vazao_massica_total = sum(vazao_massica.values())

massa_de_catalisador = (0, 25)
condicoes_iniciais = list(vazao_molar.values()) + [503, 60e5, 0]

def sistema_edos(w, y, componentes):
    # y representa as variáveis dependentes (vazão molar, por exemplo)
    # y é um vetor com as quantidades de cada componente no sistema
    condicoes_iniciais = {
        "CO": y[0],
        "H2": y[1],
        "CH3OH": y[2],
        "CH4": y[3],
        "H2O": y[4],
        "CO2": y[5],
        "N2": y[6],
        "T": y[7],
        "P": y[8],
        "Q": y[9]
    }
    
    # Chama a função de cinética que calcula as taxas de geração
    taxa_de_geracao = cinetica(condicoes_iniciais, componentes)

    # Converte o dicionário de taxas de geração em uma lista de derivadas
    dy_dt = [
        taxa_de_geracao["CO"],
        taxa_de_geracao["H2"],
        taxa_de_geracao["CH3OH"],
        taxa_de_geracao["CH4"],
        taxa_de_geracao["H2O"],
        taxa_de_geracao["CO2"],
        taxa_de_geracao["N2"],
        taxa_de_geracao["T"],
        taxa_de_geracao["P"],
        taxa_de_geracao["Q"]
    ]
    
    return dy_dt

# Resolvendo a EDO
solucao = solve_ivp(
    sistema_edos,
    massa_de_catalisador,
    condicoes_iniciais,
    args=(componentes,),
    method='RK45'  # Método Runge-Kutta de 4ª e 5ª ordem
)

# Extraindo as massas de catalisador
massa_cat = solucao.t

# ==========================
#   GRÁFICO DE COMPOSIÇÃO
# ==========================

# Extraindo as soluções para cada componente
vazao_CO = solucao.y[0]
vazao_H2 = solucao.y[1]
vazao_CH3OH = solucao.y[2]
vazao_CH4 = solucao.y[3]
vazao_H2O = solucao.y[4]
vazao_CO2 = solucao.y[5]
vazao_N2 = solucao.y[6]

# Plotando os resultados
plt.figure(figsize=(10, 6))

plt.plot(massa_cat, vazao_CO, label="CO")
plt.plot(massa_cat, vazao_H2, label="H₂")
plt.plot(massa_cat, vazao_CH3OH, label="CH₃OH")
plt.plot(massa_cat, vazao_CH4, label="CH₄")
plt.plot(massa_cat, vazao_H2O, label="H₂O")
plt.plot(massa_cat, vazao_CO2, label="CO₂")
# plt.plot(massa_cat, vazao_N2, label="N₂") - mais polui que ajuda

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Vazão Molar (mol/s)', fontsize=12)
plt.title('Soluções da EDO para cada componente', fontsize=14)
plt.legend()
plt.grid(True)

# ==========================
#   PERFIL DE PRESSÃO
# ==========================
pressao = solucao.y[8]
# Plotando os resultados
plt.figure(figsize=(10, 6))
plt.plot(massa_cat, pressao, label="P [Pa]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Pressão [Pa]', fontsize=12)
plt.title('Perfil de Pressão', fontsize=14)
plt.legend()
plt.grid(True)

# ==========================
#   PERFIL DE TEMPERATURA
# ==========================

temperatura = solucao.y[7]
# Plotando os resultados
plt.figure(figsize=(10, 6))

plt.plot(massa_cat, temperatura, label="T [K]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (kg)', fontsize=12)
plt.ylabel('Temperatura [K]', fontsize=12)
plt.title('Perfil de Temperatura', fontsize=14)
plt.legend()
plt.grid(True)

# ==========================
#   PERFIL DE CALOR
# ==========================

calor = solucao.y[9]
# Plotando os resultados
plt.figure(figsize=(10, 6))

plt.plot(massa_cat, calor, label="Q [J/s]")

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
area_secao_transversal = pi*(0.0508/2)**2  # m² - tubo de 2 in

porosidade = 0.6 # suposta - Salles (2013)
diametro_particula = 0.00311  # m - Salles (2013)
densidade_particula_catalisador = 3950  # kg_cat/m³ - Salles (2013)

densidade_do_leito = 3950 * (1-porosidade)
volume = max(massa_cat)/densidade_do_leito
altura_reator = volume/area_secao_transversal

print(f"Altura do Reator = {altura_reator:.2f} m")

# DIAMETRO
# Usando um passe triangular

diametro_externo = 0.06 # m - confirmar, considerei tubos de 2 in e e schedule 40 
pitch = 1.25*diametro_externo # distância entre os centros dos tubos considerando arranjo triangular 
diametro_casco = pitch*((numero_de_tubos/0.315)**0.5)
print(f"Diâmetro do Reator = {diametro_casco:.2f} m")

vazao_CH3OH_final = max(vazao_CH3OH) * numero_de_tubos * 3600 * componentes["CH3OH"]["MM"]*(1/1000) 
print(f"Vazão Final de Matanol = {vazao_CH3OH_final:.2f} kg/h")

# ==========================================
# GERAÇÃO DO RELATÓRIO DE SAÍDA (JSON)
# ==========================================

# 1. Isolando as vazões iniciais (mol/s) - Índice 0 da solução
F_CO_in_val    = float(solucao.y[0][0])
F_H2_in_val    = float(solucao.y[1][0])
F_CH3OH_in_val = float(solucao.y[2][0])
F_CH4_in_val   = float(solucao.y[3][0])
F_H2O_in_val   = float(solucao.y[4][0])
F_CO2_in_val   = float(solucao.y[5][0])
F_N2_in_val    = float(solucao.y[6][0])

T_in_val       = float(solucao.y[7][0])
P_in_val       = float(solucao.y[8][0])

# 2. Isolando as vazões finais (mol/s) - Índice -1 da solução
F_CO_out_val    = float(solucao.y[0][-1])
F_H2_out_val    = float(solucao.y[1][-1])
F_CH3OH_out_val = float(solucao.y[2][-1])
F_CH4_out_val   = float(solucao.y[3][-1])
F_H2O_out_val   = float(solucao.y[4][-1])
F_CO2_out_val   = float(solucao.y[5][-1])
F_N2_out_val    = float(solucao.y[6][-1])

T_out_val       = float(solucao.y[7][-1])
P_out_val       = float(solucao.y[8][-1])

# 3. Cálculos de Performance do Reator
conversao_CO = ((F_CO_in_val - F_CO_out_val) / F_CO_in_val) * 100.0 if F_CO_in_val > 0 else 0.0
conversao_CO2 = ((F_CO2_in_val - F_CO2_out_val) / F_CO2_in_val) * 100.0 if F_CO2_in_val > 0 else 0.0

# Produção de Metanol: mol/s -> kg/h total da planta
producao_metanol_tubo = F_CH3OH_out_val * 3600.0 * (componentes["CH3OH"]["MM"] / 1000.0) # kg/h por tubo
producao_metanol_total = producao_metanol_tubo * numero_de_tubos        # kg/h global

# 4. Estruturação do Dicionário de Relatório
relatorio_metanol = {
    "Parametros_Reator_Global": {
        "Numero_de_Tubos": int(numero_de_tubos),
        "Massa_Catalisador_Por_Tubo_kg": float(solucao.t[-1]),
        "Altura_Estimada_m": float(altura_reator),
        "Diametro_Casco_Estimado_m": float(diametro_casco)
    },
    "Entrada_Reator_Metanol": {
        "Temperatura_K": T_in_val,
        "Pressao_Pa": P_in_val,
        "Vazoes_Molares_Por_Tubo_mol_s": {
            "CO": F_CO_in_val,
            "H2": F_H2_in_val,
            "CH3OH": F_CH3OH_in_val,
            "CH4": F_CH4_in_val,
            "H2O": F_H2O_in_val,
            "CO2": F_CO2_in_val,
            "N2": F_N2_in_val
        }
    },
    "Saida_Reator_Metanol": {
        "Temperatura_K": T_out_val,
        "Pressao_Pa": P_out_val,
        "Vazoes_Molares_Por_Tubo_mol_s": {
            "CO": F_CO_out_val,
            "H2": F_H2_out_val,
            "CH3OH": F_CH3OH_out_val,
            "CH4": F_CH4_out_val,
            "H2O": F_H2O_out_val,
            "CO2": F_CO2_out_val,
            "N2": F_N2_out_val
        },
        "Performance_Global": {
            "Conversao_CO_percent": float(conversao_CO),
            "Conversao_CO2_percent": float(conversao_CO2),
            "Producao_Metanol_kg_h": float(producao_metanol_total)
        }
    }
}

# 4. Salvando o arquivo JSON na raiz do projeto
# (Garante que o arquivo seja salvo no mesmo diretório do script main.py)
caminho_arquivo = os.path.join(os.path.dirname(__file__), 'saida_metanol.json')

with open(caminho_arquivo, 'w', encoding='utf-8') as f:
    json.dump(relatorio_metanol, f, indent=4)

print(f"Relatório de síntese de metanol gerado com sucesso em: {caminho_arquivo}")