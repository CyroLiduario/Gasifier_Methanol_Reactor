import json
from thermodynamic import SRK, fugacidade, capacidade_calorífica, viscosidade
from math import exp
from scipy.integrate import solve_ivp

with open('./materiais.json', 'r') as materiais:
    componentes = json.load(materiais)


def cinetica(vazao_molar, temperatura, pressao, componentes, parametros):
    
    vazao_molar_total = sum(vazao_molar.values())
    vazao_massica = {}
    for componente in componentes:
        vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

    vazao_massica_total = sum(vazao_massica.values())
    
    vazao_volumetrica = vazao_molar_total * parametros["Z"]*8.314*temperatura*pressao

    rt=8.314*temperatura
    ka=1.88e8*exp(-113711/rt)
    kb=1.16e10*exp(-126523/rt)
    kc=7.08e4*exp(-68252/rt)
    KH2O=3.8e-10*exp(80876/rt)
    Kpa=exp((1.183e4/temperatura)-26.061)
    Kpb=exp((-4.773e3/temperatura)+4.672)
    Kpc=Kpa*Kpb;
    KCO=7.99e-7*exp(58700/rt)
    KCO2=1.02e-7*exp(68000/rt)
    RKH2=KH2O/(5.64e-11*exp(105600/rt))

    phi=fugacidade(vazao_massica, vazao_massica_total, parametros, componentes)
    pressoes_parciais = {}

    for componente in componentes:
        fracao_molar = vazao_molar[componente] / vazao_molar_total
        pressoes_parciais[componente] = fracao_molar * phi[componente] * (pressao / 101325)

    r1=(ka*KCO*((pressoes_parciais["CO"]*(pressoes_parciais["H2"]**1.5))-(pressoes_parciais["CH3OH"]/(Kpa*(pressoes_parciais["H2"]**0.5)))))/((1+(KCO*pressoes_parciais["CO"]))*(1+(RKH2*(pressoes_parciais["H2"]**0.5))+(KH2O*pressoes_parciais["H2O"])))
    r2=(kb*KCO2*((pressoes_parciais["CO"]*pressoes_parciais["H2"])-((pressoes_parciais["CO"]*pressoes_parciais["H2O"])/Kpb)))/((1+(KCO2*pressoes_parciais["CO2"]))*(1+(RKH2*(pressoes_parciais["H2"]**0.5))+(KH2O*pressoes_parciais["H2O"])))
    r3=(kc*KCO2*((pressoes_parciais["CO2"]*(pressoes_parciais["H2"]**1.5))-((pressoes_parciais["H2O"]*pressoes_parciais["CH3OH"])/(Kpc*(pressoes_parciais["H2"]**1.5)))))/((1+(KCO2*pressoes_parciais["CO2"]))*(1+(RKH2*(pressoes_parciais["H2"]**0.5))+(KH2O*pressoes_parciais["H2O"])))

    taxa_de_geração = {}

    taxa_de_geração["CO"] = (-r1)+(r2)
    taxa_de_geração["H2"]=(-2*r1)+(-r2)+(-3*r3)
    taxa_de_geração["CH3OH"]=(r1+r3)
    taxa_de_geração["CH4"]=(0)
    taxa_de_geração["H2O"]=(r2+r3)
    taxa_de_geração["CO2"]=(-r2)+(-r3)
    taxa_de_geração["N2"]=0

    return taxa_de_geração

def sistema_edos(w, y, temperatura, pressao, componentes, parametros):
    # y representa as variáveis dependentes (vazão molar, por exemplo)
    # y é um vetor com as quantidades de cada componente no sistema
    vazao_molar = {
        "CO": y[0],
        "H2": y[1],
        "CH3OH": y[2],
        "CH4": y[3],
        "H2O": y[4],
        "CO2": y[5],
        "N2": y[6]
    }
    
    # Chama a função de cinética que calcula as taxas de geração
    taxa_de_geracao = cinetica(vazao_molar, temperatura, pressao, componentes, parametros)

    # Converte o dicionário de taxas de geração em uma lista de derivadas
    dy_dt = [
        taxa_de_geracao["CO"],
        taxa_de_geracao["H2"],
        taxa_de_geracao["CH3OH"],
        taxa_de_geracao["CH4"],
        taxa_de_geracao["H2O"],
        taxa_de_geracao["CO2"],
        taxa_de_geracao["N2"]
    ]
    
    return dy_dt

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

alfa, aci, ai, bi, aij, a_mix, b_mix, A, B, Z = SRK(vazao_massica, vazao_massica_total, 298, 1e5,componentes)
parametros = {
    "alfa" : alfa,
    "aci" : aci,
    "ai" : ai,
    "bi" : bi,
    "aij" : aij,
    "a_mix" : a_mix,
    "b_mix" : b_mix,
    "A" : A,
    "B" : B,
    "Z" : Z
}

massa_de_catalisador = (0, 20000)
vazao_molar_inicial = list(vazao_molar.values())
# Resolvendo a EDO
solucao = solve_ivp(
    sistema_edos,
    massa_de_catalisador,
    vazao_molar_inicial,
    args=(503, 60e5, componentes, parametros),
    method='RK45'  # Método Runge-Kutta de 4ª e 5ª ordem
)

import matplotlib.pyplot as plt

# Extraindo os tempos (ou massas de catalisador)
tempo = solucao.t

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

# Mostrar gráfico
plt.show()
