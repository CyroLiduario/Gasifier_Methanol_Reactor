from math import exp
from thermodynamic import coef_fugacidade, viscosidade, SRK, capacidade_calorífica

def cinetica(condicoes_iniciais, componentes):
    
    vazao_molar = {componente: valor for componente, valor in condicoes_iniciais.items() if componente not in ("P", "T")}
    temperatura = condicoes_iniciais["T"]
    pressao = condicoes_iniciais["P"]

    print(f'Pressão: {pressao} Pa')
    print(f'Temperatura: {temperatura} K')

    vazao_molar_total = sum(vazao_molar.values())
    vazao_massica = {}
    for componente in componentes:
        vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

    vazao_massica_total = sum(vazao_massica.values())
    
    print(f"Vazão Mássica: {vazao_massica_total}")
    
    alfa, aci, ai, bi, aij, a_mix, b_mix, A, B, Z = SRK(vazao_molar, vazao_molar_total, temperatura, pressao, componentes)
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

    vazao_volumetrica = vazao_molar_total * parametros["Z"]*8.314*temperatura/pressao

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

    phi=coef_fugacidade(vazao_molar, vazao_molar_total, parametros, componentes)
    pressoes_parciais = {}

    for componente in componentes:
        fracao_molar = vazao_molar[componente] / vazao_molar_total
        pressoes_parciais[componente] = fracao_molar * phi[componente] * (pressao / 101325)
        #pressoes_parciais[componente] = fracao_molar * phi * (pressao / 101325)

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

    taxa_de_geração["P"] = perda_de_carga(componentes, vazao_massica, vazao_massica_total, temperatura, vazao_volumetrica)
    taxa_de_geração["T"] = balanco_energia(componentes, temperatura, vazao_molar, r1, r2, r3)

    return taxa_de_geração

def balanco_energia(componentes, temperatura, vazao_molar, r1, r2, r3):
    # ==============================
    #   BALANÇO DE ENERGIA
    # ==============================

    # entalpias padrão de reação
    dHp_r1 = (-90130) # reação de hidrogenação do monóxido de carbono
    dHp_r2 = 41470 # water gas shift
    dHp_r3 = (-48970) # hidrogenação do dióxido

    #   capacidades caloríficas por reação
    cp = capacidade_calorífica(componentes, temperatura)
    dCp_r1 = (-cp["CO"])-(2*(cp["H2"]))+cp["CH3OH"]
    dCp_r2 = (-cp["CO2"])-cp["H2"]+cp["CO"]+cp["H2O"]
    dCp_r3 = (-cp["CO2"])-(3*cp["H2"])+cp["CH3OH"]+cp["H2O"]

    # entalpia de reação na temperatura de operação
    dH_r1 = dHp_r1 + (dCp_r1*(temperatura-298))
    dH_r2 = dHp_r2 + (dCp_r2*(temperatura-298))
    dH_r3 = dHp_r3 + (dCp_r3*(temperatura-298))

    # calor específico da corrente
    sum_FiCpi = sum(vazao_molar[k] * cp[k] for k in vazao_molar if k not in ("T", "P"))

    sum_dH = (r1*(-dH_r1))+(r2*(-dH_r2))+(r3*(-dH_r3))
    Q= 0
    #dTdW = (Q+sum_dH)/sum_FiCpi
    dTdW = 0
    return dTdW

def perda_de_carga(componentes, vazao_massica, vazao_massica_total, temperatura, vazao_volumetrica):
    # ================================
    #   PERDA DE CARGA POR ERGUN
    # ================================

    # === TO-DO: finalizar código de perda de carga ===
    # densidade_catalizador = 3950 #kg/m³
    # area_transversal = 0.002 #m²
    # visc_mistura = viscosidade(componentes, vazao_massica, vazao_massica_total, temperatura)
    # densidade_mistura = vazao_massica_total/vazao_volumetrica
    # velocidade = vazao_volumetrica/area_transversal #m²
    # porosidade = 0.55 #suposto
    # diametro_particula = 0.00311 #diametro do catalizador (artigo)
    # taxa_massica = densidade_mistura*velocidade
    # reynolds = diametro_particula*densidade_mistura*velocidade/(visc_mistura*(1-porosidade))
    # dPdW = 0

    # ===============================

    area_secao_transversal = 0.002  # m² (ajuste caso este valor represente o diâmetro)
    porosidade = 0.55
    diametro_particula = 0.00311  # m
    densidade_particula_catalisador = 3950  # kg_cat/m³
    
    # 1. Ajuste de Unidades para o SI (garantindo que massa esteja em kg para a Ergun)
    # Se o MM no json for g/mol, vazao_massica está em g/s. Convertendo para kg/s:
    vazao_massica_kg_s = vazao_massica_total
    
    # Viscosidade da mistura (garanta que a função retorne em Pa.s ou kg/(m.s))
    visc_mistura = viscosidade(componentes, vazao_massica, vazao_massica_total, temperatura)
    
    # 2. Propriedades Fluidodinâmicas
    velocidade = vazao_volumetrica / area_secao_transversal # m/s
    densidade_gas = vazao_massica_kg_s / vazao_volumetrica # kg/m³
    
    # 3. Equação de Ergun (Calcula dP/dz em Pa/m)
    termo_viscoso = (150 * visc_mistura * (1 - porosidade)**2) / ((diametro_particula**2) * (porosidade**3))
    termo_inercial = (1.75 * densidade_gas * (1 - porosidade)) / (diametro_particula * (porosidade**3))
    
    dP_dz = - (termo_viscoso * velocidade + termo_inercial * (velocidade**2))
    
    # 4. Conversão de domínio (z para W)
    # dW_dz = Massa de catalisador por metro de tubo (kg_cat / m)
    dW_dz_kg = area_secao_transversal * densidade_particula_catalisador * (1 - porosidade)
    
    # Como o nosso solver está integrando em GRAMAS de catalisador, convertemos:
    dW_dz_g = dW_dz_kg * 1000 # g_cat / m
    
    # Queda de pressão por grama de catalisador (Pa / g_cat)
    dPdW = dP_dz / dW_dz_g
    print(dPdW)
    return dPdW