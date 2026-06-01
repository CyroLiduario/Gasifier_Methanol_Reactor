from math import exp
from thermodynamic import coef_fugacidade, viscosidade

def cinetica(condicoes_iniciais, componentes, parametros):
    
    vazao_molar = {componente: valor for componente, valor in condicoes_iniciais.items() if componente not in ("P", "T")}
    temperatura = condicoes_iniciais["T"]
    pressao = condicoes_iniciais["P"]

    vazao_molar_total = sum(vazao_molar.values())
    vazao_massica = {}
    for componente in componentes:
        vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

    vazao_massica_total = sum(vazao_massica.values())
    
    vazao_volumetrica = vazao_molar_total * parametros["Z"]*8.314*temperatura*pressao

    # PERDA DE CARGA POR ERGUN
    # === TO-DO: finalizar código de perda de carga ===
    densidade_catalizador = 3950
    visc_mistura = viscosidade(componentes, vazao_massica, vazao_massica_total, temperatura)
    densidade_mistura = vazao_massica_total/vazao_volumetrica
    velocidade = vazao_volumetrica*(0.002) #m²
    porosidade = 0.55 #suposto
    diametro_particula = 0.00311 #diametro do catalizador (artigo)
    taxa_massica = densidade_mistura*velocidade
    reynolds = diametro_particula*densidade_mistura*velocidade/(visc_mistura*(1-porosidade))

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

    phi=coef_fugacidade(vazao_massica, vazao_massica_total, parametros, componentes)
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

    taxa_de_geração["P"] = 0
    taxa_de_geração["T"] = 0

    return taxa_de_geração

