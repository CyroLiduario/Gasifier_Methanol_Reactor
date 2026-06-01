import json
import matplotlib.pyplot as plt
from kinetic_model import cinetica
from scipy.integrate import solve_ivp

with open('./materiais.json', 'r') as materiais:
    componentes = json.load(materiais)

# Confirmar se as vazões estão em kmol/h
# Arbitrar um número de tubos (dividir a vazão total pelo número de tubos)
vazao_molar = {
    "CO" : 945,
    "H2" : 1850,
    "CH3OH" : 18,
    "CH4" : 29,
    "H2O" : 26,
    "CO2" : 483,
    "N2" : 846
}

vazao_massica = {}
for componente in componentes:
    vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

vazao_massica_total = sum(vazao_massica.values())

massa_de_catalisador = (0, 20000)
condicoes_iniciais = list(vazao_molar.values()) + [503, 60e5]

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
        "P": y[8]
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
        taxa_de_geracao["P"]
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

# Extraindo os tempos (ou massas de catalisador)
tempo = solucao.t

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

plt.plot(tempo, vazao_CO, label="CO")
plt.plot(tempo, vazao_H2, label="H₂")
plt.plot(tempo, vazao_CH3OH, label="CH₃OH")
plt.plot(tempo, vazao_CH4, label="CH₄")
plt.plot(tempo, vazao_H2O, label="H₂O")
plt.plot(tempo, vazao_CO2, label="CO₂")
plt.plot(tempo, vazao_N2, label="N₂")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (g)', fontsize=12)
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

plt.plot(tempo, pressao, label="P [Pa]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (g)', fontsize=12)
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

plt.plot(tempo, temperatura, label="T [K]")

# Personalização do gráfico
plt.xlabel('Massa de Catalisador (g)', fontsize=12)
plt.ylabel('Temperatura [K]', fontsize=12)
plt.title('Perfil de Temperatura', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()