# Define thermodynamic equations that will be used

import numpy as np

def SRK(massa, massa_total, temperatura, pressão, componentes):
        
    y = {}
    for componente in componentes:
        y[componente] = massa[componente]/massa_total
    
    R = 8.314 # J/mol.K

    alfa = {}
    aci = {}
    ai = {}
    bi = {}

    for componente in componentes:
        alfa[componente]=(1+(0.480+(1.574*componentes[componente]["w"])-(0.176*componentes[componente]["w"]**2))*(1-((temperatura/componentes[componente]["Tc"])**0.5)))**2
        aci[componente]=0.42748*(((R**2)*(componentes[componente]["Tc"]**2.5))/componentes[componente]["Pc"])
        ai[componente]=alfa[componente]*aci[componente]
        bi[componente]=0.008664*R*componentes[componente]["Tc"]/componentes[componente]["Pc"]
    
    #Regra de mistura
    aij = {}
    a_mix = 0
    b_mix = 0
    for componente1 in componentes:
        for componente2 in componentes:
            aij[(componente1, componente2)] = y[componente1]*y[componente2]*((ai[componente1]*ai[componente2])**0.5) #* (1 - kij.get((i, j), kij.get((j, i), 0))) - parametro considerado 0 por simplificação
            a_mix += aij[(componente1, componente2)] 
        b_mix += y[componente1]*bi[componente1]

    
    A=(a_mix*pressão)/((R**2)*(temperatura**2))
    B=(b_mix*pressão)/(R*temperatura)
    
    #coeffs = [1, -(1 - B), A - 2*B - 3*B**2, -(A*B - B**2 - B**3)]
    coeffs = [1, -1, (A - B - B**2), -(A * B)]
    Z_roots = np.roots(coeffs)
    Z = np.real(Z_roots[np.isreal(Z_roots)])[0]

    return alfa, aci, ai, bi, aij, a_mix, b_mix, A, B, Z

def fugacidade(massa, massa_total, parametros, componentes):
    y = {}
    for componente in componentes:
        y[componente] = massa[componente]/massa_total
    
    phi = {}
    for componente1 in componentes:
        term1 = parametros["bi"][componente1] * (parametros["Z"] - 1) / parametros["b_mix"]
        term2 = np.log(parametros["Z"] - parametros["B"])
        term3 = parametros["A"] / parametros["B"]
        sum_term4 = 0
        for componente2 in componentes:
            sum_term4 += y[componente2]*np.sqrt(parametros["ai"][componente1]*parametros["ai"][componente2])
        
        term4 = (2* sum_term4)/parametros["a_mix"]
        term5 = np.log((parametros["Z"]+parametros["B"])/parametros["Z"])
        ln_phi_i = term1 - term2 - term3 * (term4 - (parametros["bi"][componente1]/parametros["b_mix"]))*term5
        phi[componente1] = np.exp(ln_phi_i)
    
    return phi

def capacidade_calorífica(componentes, temperatura):
    R = 8.314
    cp = {}
    for componente in componentes:
        cp[componente] = R*(componente["A"]+(componente["B"]*temperatura)+(componente["C"]*(temperatura**2))+(componente["D"]*(temperatura**-2))+(componente["E"]*(temperatura**3)))
    
    return cp

def viscosidade(componentes, massa, massa_total,temperatura):
    y = {}
    visc = {}
    visc_mix = 0

    for componente in componentes:
        y[componente] = massa[componente]/massa_total
        visc[componente] = (componentes[componente]["constantes_de_viscosidade"][0]*(temperatura**componentes[componente]["constantes_de_viscosidade"][1]))/(1+(componentes[componente]["constantes_de_viscosidade"][2]/temperatura)+(componentes[componente]["constantes_de_viscosidade"][3]/temperatura**2))
        visc_mix += y[componente]*visc[componente]

    return visc_mix