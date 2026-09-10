from math import exp, pi
from thermodynamic import capacidade_calorífica, viscosidade

def cinetica_reformer(condicoes_iniciais, componentes, diametro_reator):
    # Baseado nos estudos de Ramos et al. (2024) e García-Vargas et al. (2015)
    '''
        Catalyst: Ni-Mg/β-SiC
        Dry Reforming: CH4 + CO2 ↔ 2 CO + 2 H2
        Steam Reforming: CH4 + H2O ↔ CO+ 3 H2
        Water-gas shift: CO + H2O ↔ CO2 + H2
        Reverse Methanation: CH4 + 2 H2O ↔ CO2 + 4 H2 - desconsiderada no estudo de Ramos et al. (2024)
    '''

    vazao_molar = {componente: valor for componente, valor in condicoes_iniciais.items() if componente not in ("P", "T", "Q")}
    temperatura = condicoes_iniciais["T"]
    pressao = condicoes_iniciais["P"]
    vazao_molar_total = sum(vazao_molar.values())

    # vazao massica
    vazao_massica = {}
    for componente in componentes:
        vazao_massica[componente] = componentes[componente]["MM"] * vazao_molar[componente]

    vazao_massica_total = sum(vazao_massica.values())
    vazao_volumetrica = vazao_molar_total *8.314*temperatura/pressao #considerando gás ideal pela baixa pressão

    pressoes_parciais = {}

    for componente in componentes:
        fracao_molar = vazao_molar[componente] / vazao_molar_total
        pressoes_parciais[componente] = fracao_molar * (pressao / 1000) # kPa
    
    #Rg = 8.314 # J/(mol.K)
    Rg = 0.008314 # kJ/mol.K

    # Dry Reforming: CH4 + CO2 ↔ 2 CO + 2 H2
    k_dr = 70.99*exp(-77.820/(Rg*temperatura)) 
    K_dr = 6.780e18*exp(-31230/temperatura)
    r_dr = k_dr*pressoes_parciais["CH4"]*(1-(((pressoes_parciais["CO"]**2)*(pressoes_parciais["H2"]**2))/(pressoes_parciais["CH4"]*pressoes_parciais["CO2"]*K_dr)))

    # Steam Reforming: CH4 + H2O ↔ CO+ 3 H2
    k_sr = 85.77*exp(-74.720/(Rg*temperatura)) 
    K_sr = 1.198e17*exp(-26830/temperatura)
    r_sr = k_sr*pressoes_parciais["CH4"]*(1-((pressoes_parciais["CO"]*(pressoes_parciais["H2"]**3))/(pressoes_parciais["CH4"]*pressoes_parciais["H2O"]*K_sr)))

    # Water-gas shift: CO + H2O ↔ CO2 + H2
    k_wgs = 149.92*exp(-54.260/(Rg*temperatura)) 
    K_wgs = 10**((2078/temperatura)-2.029)
    r_wgs = k_wgs*(((pressoes_parciais["CO"]*pressoes_parciais["H2O"])/pressoes_parciais["H2"])-(pressoes_parciais["CO2"]/K_wgs))
    
    # 5. Balanço Estequiométrico das Espécies (mol / kg_cat . s)
    taxa_de_geracao = {}
    taxa_de_geracao["H2O"] = (0)+(-r_sr)+(-r_wgs)
    taxa_de_geracao["CO"] = (2*r_dr)+(r_sr)+(-r_wgs)
    taxa_de_geracao["CO2"] = (-r_dr)+(0)+(r_wgs)
    taxa_de_geracao["H2"] = (2*r_dr)+(3*r_sr)+(r_wgs)
    taxa_de_geracao["CH4"] = (-r_dr)+(-r_sr)+(0)
    taxa_de_geracao["N2"] = 0

    taxa_de_geracao["P"] = perda_de_carga(componentes, vazao_massica, vazao_massica_total, temperatura, vazao_volumetrica, diametro_reator)
    # taxa_de_geracao["P"] = 0

    dTdW, dQdW = balanco_energia(componentes, temperatura, vazao_molar, r_dr, r_sr, r_wgs)
    
    taxa_de_geracao["T"] = dTdW
    taxa_de_geracao["Q"] = dQdW
    
    return taxa_de_geracao

def balanco_energia(componentes, temperatura, vazao_molar, r_dr, r_sr, r_wgs):
    dH_dr = 247000 #J/mol
    dH_sr = 206000 #J/mol
    dH_wgs = -41000 #J/mol

    cp = capacidade_calorífica(componentes, temperatura)

    # Dry Reforming: CH4 + CO2 ↔ 2 CO + 2 H2
    dCp_r_dr = (2*cp["CO"])+(2*(cp["H2"]))+(-cp["CH4"])+(-cp["CO2"])

    # Steam Reforming: CH4 + H2O ↔ CO+ 3 H2
    dCp_r_sr = (cp["CO"])+(3*cp["H2"])+(-cp["CH4"])+(-cp["H2O"])

    # Water-gas shift: CO + H2O ↔ CO2 + H2
    dCp_r_wgs = (cp["CO2"])+(cp["H2"])+(-cp["CO"])+(-cp["H2O"])

    # Entalpia de reação na temperatura de operação - J/mol.K
    dH_dr = dH_dr + (dCp_r_dr*(temperatura-298))
    dH_sr = dH_sr + (dCp_r_sr*(temperatura-298))
    dH_wgs = dH_wgs + (dCp_r_wgs*(temperatura-298))

    sum_FiCpi = sum(vazao_molar[k] * cp[k] for k in vazao_molar if k not in ("T", "P", "Q"))

    sum_dH = (r_dr*(-dH_dr))+(r_sr*(-dH_sr))+(r_wgs*(-dH_wgs))
    #Q = 0 # modelando um processo adiabático - sem troca térmica
    Q = -sum_dH
    dTdW = (Q+sum_dH)/sum_FiCpi
    return dTdW, Q

def perda_de_carga(componentes, vazao_massica, vazao_massica_total, temperatura, vazao_volumetrica, diametro_reator):
    # ================================
    #   PERDA DE CARGA POR ERGUN
    # ================================
    area_secao_transversal = pi*(diametro_reator/2)**2  # m² - tubo de 2 in

    porosidade = 0.6 # suposta - Salles (2013)
    diametro_particula = 0.001  # m - GARCÍA-VARGAS et al. (2015)
    densidade_particula_catalisador = (0.05*8900)+(0.008*1738)+((1-0.05-0.008)*3210)  # kg_cat/m³ - Salles (2013)
    
    # Convertendo a vazão mássica para kg/s:
    vazao_massica_kg_s = vazao_massica_total/1000
    
    # Viscosidade da mistura (garanta que a função retorne em Pa.s ou kg/(m.s))
    visc_mistura = viscosidade(componentes, vazao_massica, vazao_massica_total, temperatura)

    velocidade = vazao_volumetrica / area_secao_transversal # m/s
    
    densidade_gas = vazao_massica_kg_s / vazao_volumetrica # kg/m³
    
    # 3. Equação de Ergun (Calcula dP/dz em Pa/m)
    termo_viscoso = (150 * visc_mistura * (1 - porosidade)**2) / ((diametro_particula**2) * (porosidade**3))
    termo_inercial = (1.75 * densidade_gas * (1 - porosidade)) / (diametro_particula * (porosidade**3))
    
    dP_dz = - (termo_viscoso * velocidade + termo_inercial * (velocidade**2))
    
    # Conversão de domínio (z para W)
    # dW_dz = Massa de catalisador por metro de tubo (kg_cat / m)
    dW_dz_kg = area_secao_transversal * densidade_particula_catalisador * (1 - porosidade)
        
    # Queda de pressão por grama de catalisador (Pa / kg_cat)
    dPdW = dP_dz / dW_dz_kg

    return dPdW