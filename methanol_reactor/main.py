import json
import matplotlib.pyplot as plt
from math import pi
from kinetic_model import cinetica
from scipy.integrate import solve_ivp

with open('./materiais.json', 'r') as materiais:
    componentes = json.load(materiais)

# Confirmar se as vazões estão em kmol/h
# Arbitrar um número de tubos (dividir a vazão total pelo número de tubos)

vazao_molar = { #mol/s
    "CO" : 138.0,
    "H2" : 114.3,
    "CH3OH" : 0,
    "CH4" : 4.5,
    "H2O" : 0,
    "CO2" : 95.0,
    "N2" : 351.4
}

vazao_molar = {component:vazao*(3600/1000) for component,vazao in vazao_molar.items()}

""" vazao_molar = { #kmol/h
    "CO" : 3418.4,
    "H2" : 6788.6,
    "CH3OH" : 65,
    "CH4" : 303.1,
    "H2O" : 42.8,
    "CO2" : 42.8,
    "N2" : 3048.7
} """
#print(vazao_molar)

# transformando em mol/s

vazao_molar = { componente: vazao_molar_comp*(1000/3600) for componente, vazao_molar_comp in vazao_molar.items()}
#print(vazao_molar)

# vazao por tubo
numero_de_tubos = 1000
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
plt.plot(massa_cat, vazao_N2, label="N₂")

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

#print(vazao_CH3OH)
vazao_CH3OH_final = max(vazao_CH3OH) * numero_de_tubos * 3600 * componentes["CH3OH"]["MM"]*(1/1000) 
print(f"Vazão Final de Matanol = {vazao_CH3OH_final:.2f} kg/h")