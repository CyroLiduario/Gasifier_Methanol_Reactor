import os
import json
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import cinetica_reformer
from scipy.optimize import minimize

# ==========================================
# 1. LEITURA DOS DADOS DE SAÍDA DO GASEIFICADOR
# ==========================================
quantidade_reatores = 1
diametro_reator = 2.5
caminho_saida_gas = os.path.join(os.path.dirname(__file__), '..', 'gasifier', 'saida_gaseificador.json')

with open(caminho_saida_gas, 'r', encoding='utf-8') as f:
    dados_in = json.load(f)

saida = dados_in["Saida_Gaseificador"]
vazoes_in = dict(saida["Vazoes_Molares_mol_s"])  

# Removendo O2 e Tar (são zerados na saída do gaseificador)
vazoes_in.pop("O2", None)
vazoes_in.pop("Tar", None)

T_in = saida["Temperatura_K"]
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
vazao_vapor_adicional = 600  # mol/s
vazoes_in["H2O"] = vazoes_in["H2O"] + vazao_vapor_adicional

# Injeção de CH4 adicional
vazao_CH4_adicional = 600 # mol/s
vazoes_in["CH4"] = vazoes_in["CH4"] + vazao_CH4_adicional

vazoes_in = {componente: vazao_molar_comp/quantidade_reatores for componente, vazao_molar_comp in vazoes_in.items()}
#numero_de_tubos = 1
#vazoes_in = {componente: vazao_molar_comp/numero_de_tubos for componente, vazao_molar_comp in vazoes_in.items()}

Massa_Catalisador_Total = (0, 1000) # kg de catalisador

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
        "T": y[6],
        "P": y[7],
        "Q": y[8]
    }
    
    # Chama a função de cinética que calcula as taxas de geração
    taxa_de_geracao = cinetica_reformer(condicoes_iniciais, componentes, diametro_reator)

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
def avaliar_reformer(x, componentes, diametro_reator, Massa_Catalisador_Total, quantidade_reatores, vazoes_brutas, P_in):
    # Variáveis de decisão testadas pelo algoritmo
    vazao_vapor_adicional = x[0]
    vazao_CH4_adicional = x[1]
    T_in_opt = x[2]
    
    # Reconstrói as vazões de entrada (dividindo pelo número de reatores paralelos)
    vazoes_in_opt = dict(vazoes_brutas)
    vazoes_in_opt.pop("O2", None)
    vazoes_in_opt.pop("Tar", None)
    
    vazoes_in_opt["H2O"] += vazao_vapor_adicional
    vazoes_in_opt["CH4"] += vazao_CH4_adicional
    vazoes_in_opt = {comp: v / quantidade_reatores for comp, v in vazoes_in_opt.items()}
    
    y0_opt = [
        vazoes_in_opt.get("CO", 0), 
        vazoes_in_opt.get("H2", 0), 
        vazoes_in_opt.get("CH4", 0), 
        vazoes_in_opt.get("H2O", 0), 
        vazoes_in_opt.get("CO2", 0), 
        vazoes_in_opt.get("N2", 0), 
        T_in_opt, 
        P_in, 
        0.0 # Q
    ]
    
    # Integração rápida usando BDF (Ideal para algoritmos de otimização)
    try:
        sol = solve_ivp(
            sistema_edos,
            Massa_Catalisador_Total,
            y0_opt,
            args=(componentes, diametro_reator),
            method='BDF', 
            rtol=1e-3, atol=1e-6
        )
        
        # Extrai os valores finais
        F_CO_out = sol.y[0][-1]
        F_H2_out = sol.y[1][-1]
        P_out = sol.y[7][-1]
        
        razao_out = F_H2_out / F_CO_out if F_CO_out > 0 else 0
        
    except:
        # Se o integrador falhar por condições termodinâmicas extremas, retorna uma penalidade gigante
        return 1e6

    # ----------------------------------------------------
    # FUNÇÕES OBJETIVO
    # ----------------------------------------------------
    # F1: Erro quadrático da razão H2/CO alvo (Peso Gigante para forçar a estequiometria)
    F1 = (razao_out - 2.0)**2
    
    # F2: Custo de injeção de insumos (Escalado para não ofuscar o F1)
    F2 = (vazao_vapor_adicional + vazao_CH4_adicional) / 1000.0
    
    # Restrição: Penalidade por queda de pressão (se a pressão cair mais de 0.5 bar)
    queda_pressao_bar = (P_in - P_out) / 100000.0
    penalidade_pressao = 5000.0 if queda_pressao_bar > 0.5 else 0.0
    
    # Custo Total Ponderado
    custo_total = (1000.0 * F1) + (1.0 * F2) + penalidade_pressao
    
    return custo_total

# Limites Operacionais (Bounds)
# (Vapor: 0 a 1000 mol/s) | (Metano: 0 a 1000 mol/s) | (T_in: 800 a 1100 °C convertidos para K)
bounds_moo = [
    (0.0, 1000.0), 
    (0.0, 1000.0), 
    (1138.9, 1373.15)
]

# Chute inicial razoável para o algoritmo começar a procurar
x0 = [300.0, 300.0, 1138.0]

print("Iniciando Otimização Multiobjetivo (Isso pode levar alguns minutos)...")

resultado_opt = minimize(
    avaliar_reformer, 
    x0, 
    args=(componentes, diametro_reator, Massa_Catalisador_Total, quantidade_reatores, dados_in["Saida_Gaseificador"]["Vazoes_Molares_mol_s"], P_in),
    method='SLSQP', 
    bounds=bounds_moo, 
    options={'disp': True, 'maxiter': 50, 'eps': 1.0}
)

if resultado_opt.success:
    vapor_opt, metano_opt, T_in_opt = resultado_opt.x
    print("\n==========================================")
    print("        PONTO ÓTIMO ENCONTRADO")
    print("==========================================")
    print(f"Vapor Adicional: {vapor_opt:.2f} mol/s")
    print(f"Metano Adicional: {metano_opt:.2f} mol/s")
    print(f"Temperatura de Entrada: {T_in_opt:.2f} K ({T_in_opt - 273.15:.2f} °C)")
    
    vazao_vapor_adicional = vapor_opt
    vazao_CH4_adicional = metano_opt
    T_in = T_in_opt

    # --- RECONSTRUÇÃO DA SOLUÇÃO ÓTIMA GLOBAL ---
    vazoes_in_opt = dict(dados_in["Saida_Gaseificador"]["Vazoes_Molares_mol_s"])
    vazoes_in_opt.pop("O2", None)
    vazoes_in_opt.pop("Tar", None)
    vazoes_in_opt["H2O"] += vapor_opt
    vazoes_in_opt["CH4"] += metano_opt
    
    # Atualiza vazão da linha pro relatório e pro reator paralelo
    vazoes_in = vazoes_in_opt.copy() 
    vazoes_in_opt = {comp: v / quantidade_reatores for comp, v in vazoes_in_opt.items()}
    
    y0_opt = [
        vazoes_in_opt.get("CO", 0), vazoes_in_opt.get("H2", 0), vazoes_in_opt.get("CH4", 0), 
        vazoes_in_opt.get("H2O", 0), vazoes_in_opt.get("CO2", 0), vazoes_in_opt.get("N2", 0), 
        T_in_opt, P_in, 0.0
    ]
    
    # Roda a simulação definitiva para extrair os gráficos e o JSON
    solucao = solve_ivp(
        sistema_edos, Massa_Catalisador_Total, y0_opt,
        args=(componentes, diametro_reator), method='BDF', rtol=1e-4, atol=1e-6
    )
else:
    print("A otimização não convergiu. Verifique os limites ou o chute inicial.")

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
        "Temperatura_Entrada_K": float(T_in),
        "Vazoes_Molares": vazoes_in
    },
    "Saida": {
        "Vazoes_Molares": {
            "CO": float(F_CO_out*quantidade_reatores), 
            "H2": float(F_H2_out*quantidade_reatores), 
            "CH4": float(F_CH4_out*quantidade_reatores),
            "H2O": float(F_H2O_out*quantidade_reatores), 
            "CO2": float(F_CO2_out*quantidade_reatores), 
            "N2": float(F_N2_out*quantidade_reatores)
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
# 7. VISUALIZAÇÃO DOS PERFIS
# ==========================================

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
ax1.legend(loc='center right')

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
diametro_particula = 0.001  # m - GARCÍA-VARGAS et al. (2015)
densidade_particula_catalisador = (0.05*8900)+(0.008*1738)+((1-0.05-0.008)*3210)  # kg_cat/m³ - Salles (2013)

densidade_do_leito = densidade_particula_catalisador * (1-porosidade)
volume = max(massa_w)/densidade_do_leito
altura_reator = volume/area_secao_transversal

print(f"Altura do Reator = {altura_reator:.2f} m")

# DIAMETRO
# Usando um passe triangular

# diametro_externo = 0.06 # m - confirmar, considerei tubos de 2 in e e schedule 40 
# pitch = 1.25*diametro_externo # distância entre os centros dos tubos considerando arranjo triangular 

print(f"Diâmetro do Reator = {diametro_reator:.2f} m")

# ==========================================
# 8. GRÁFICO DE SUPERFÍCIE (Análise de Sensibilidade 3D)
# ==========================================
print("\nGerando gráfico de superfície (isso pode levar alguns instantes)...")

# Cria uma malha de 15x15 pontos ao redor do ponto ótimo
vapor_range = np.linspace(max(0, vapor_opt - 300), vapor_opt + 300, 15)
metano_range = np.linspace(max(0, metano_opt - 300), metano_opt + 300, 15)
V, M = np.meshgrid(vapor_range, metano_range)
Razao_Matriz = np.zeros_like(V)

# O loop avalia a razão H2/CO para cada combinação de Vapor e Metano
for i in range(V.shape[0]):
    for j in range(V.shape[1]):
        v_test = V[i, j]
        m_test = M[i, j]
        
        # Reconstrói a entrada apenas para o teste
        v_test_dict = dict(dados_in["Saida_Gaseificador"]["Vazoes_Molares_mol_s"])
        v_test_dict.pop("O2", None)
        v_test_dict.pop("Tar", None)
        v_test_dict["H2O"] += v_test
        v_test_dict["CH4"] += m_test
        v_test_dict = {comp: val / quantidade_reatores for comp, val in v_test_dict.items()}
        
        y0_test = [
            v_test_dict.get("CO", 0), v_test_dict.get("H2", 0), v_test_dict.get("CH4", 0), 
            v_test_dict.get("H2O", 0), v_test_dict.get("CO2", 0), v_test_dict.get("N2", 0), 
            T_in_opt, P_in, 0.0
        ]
        
        sol_test = solve_ivp(
            sistema_edos, Massa_Catalisador_Total, y0_test,
            args=(componentes, diametro_reator), method='BDF', rtol=1e-3, atol=1e-6
        )
        
        F_CO_test = sol_test.y[0][-1]
        F_H2_test = sol_test.y[1][-1]
        Razao_Matriz[i, j] = F_H2_test / F_CO_test if F_CO_test > 0 else 0

# Plotagem da Superfície 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(V, M, Razao_Matriz, cmap='viridis', edgecolor='none', alpha=0.9)
ax.scatter(vapor_opt, metano_opt, H2_CO_out, color='red', s=100, label='Ponto Ótimo', zorder=5)

ax.set_title(f'Superfície de Resposta: Razão H₂/CO (T = {T_in_opt - 273.15:.1f} °C)')
ax.set_xlabel('Injeção de Vapor d\'água (mol/s)')
ax.set_ylabel('Injeção de Metano (mol/s)')
ax.set_zlabel('Razão H₂/CO de Saída')
fig.colorbar(surf, shrink=0.5, aspect=0.5, label='Razão H₂/CO')
ax.legend()
plt.show()

# ==========================================
# 9. ANÁLISE DE SENSIBILIDADE: Razão H2/CO vs Temperatura
# ==========================================
print("\nGerando perfil de sensibilidade térmica 2D...")

# Definindo a faixa de temperatura (ex: 800 °C a 1100 °C)
temperaturas_C = np.linspace(700, 1140, 25)
temperaturas_K = temperaturas_C + 273.15
razoes_H2_CO = []

# Mantendo as injeções ótimas fixas (calculadas na etapa 5)
vazoes_in_base = dict(dados_in["Saida_Gaseificador"]["Vazoes_Molares_mol_s"])
vazoes_in_base.pop("O2", None)
vazoes_in_base.pop("Tar", None)
vazoes_in_base["H2O"] += vapor_opt
vazoes_in_base["CH4"] += metano_opt

# Divide a corrente para o reator paralelo
vazoes_in_base = {comp: val / quantidade_reatores for comp, val in vazoes_in_base.items()}

# Loop de simulação para cada temperatura
for T_teste in temperaturas_K:
    y0_teste = [
        vazoes_in_base.get("CO", 0), vazoes_in_base.get("H2", 0), vazoes_in_base.get("CH4", 0), 
        vazoes_in_base.get("H2O", 0), vazoes_in_base.get("CO2", 0), vazoes_in_base.get("N2", 0), 
        T_teste, P_in, 0.0
    ]
    
    # Resolve a EDO para a temperatura testada (uso do solver BDF para rapidez)
    sol_teste = solve_ivp(
        sistema_edos, Massa_Catalisador_Total, y0_teste,
        args=(componentes, diametro_reator), method='BDF', rtol=1e-3, atol=1e-6
    )
    
    F_CO_final = sol_teste.y[0][-1]
    F_H2_final = sol_teste.y[1][-1]
    
    razao = F_H2_final / F_CO_final if F_CO_final > 0 else 0
    razoes_H2_CO.append(razao)

# Plotagem do Gráfico 2D
plt.figure(figsize=(10, 6))
plt.plot(temperaturas_C, razoes_H2_CO, 'b-', linewidth=2.5, label='Curva de Equilíbrio Térmico')

# Marcações de referência
plt.axvline(x=(T_in_opt - 273.15), color='red', linestyle='--', label=f'Ponto Ótimo ({T_in_opt - 273.15:.1f} °C)')
plt.axhline(y=2.0, color='green', linestyle=':', linewidth=2, label='Meta Estequiométrica (2.0)')

plt.title('Sensibilidade da Razão H₂/CO com a Temperatura de Entrada', fontsize=14)
plt.xlabel('Temperatura de Entrada (°C)', fontsize=12)
plt.ylabel('Razão H₂/CO de Saída', fontsize=12)
plt.legend()
plt.grid(True)
plt.show()