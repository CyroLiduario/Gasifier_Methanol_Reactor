import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from kinetic_model import gasifier_odes
import json

# ==========================================
# 1. PREPARAÇÃO DO AMBIENTE E CONSTANTES
# ==========================================
with open('./gasifier/materials.json', 'r') as f:
    materiais = json.load(f)

# Parâmetros Iniciais (Idênticos ao seu main)
Rg = 8.314
P0 = 3.0 * 101325.0
T0 = 300.0
reator_params = {
    "L": 12.0, "D": 2.0, "eps": 0.5, "phi": 2.0 / 3.0, "em": 0.9,
    "d_p": 0.042544, "rhochar0": 152.5, "sigma": 5.67051e-8, "T_w": 1050.0
}

area_secao = (math.pi * reator_params["D"]**2) / 4.0
F_s = 101104.0 / 5.0
YM = 0.0872
rhom_s0 = 119.0
U_s0 = (F_s / 3600.0) / (rhom_s0 * area_secao)
rho_B0 = rhom_s0 * (1.0 - YM)
rho_char0 = 0.0
rho_M0 = rhom_s0 * YM

ER = 0.18
# (Para simplificar o snippet, os cálculos de mcomb_bag, var_bag e ER_stoic entram aqui)
# ER_stoic = ... 
# F_g = ...
rhom_g0_ar = 101325.0 / (287.058 * 300.0)
U_g0 = (51480.0 / 3600.0) / (rhom_g0_ar * area_secao) # Exemplo numérico de velocidade de ar

# ==========================================
# 2. LOOP DE ANÁLISE DE SENSIBILIDADE
# ==========================================
# Vamos testar a fração de vapor de 0% (0.0) até 50% (0.5) em 11 passos (0%, 5%, 10%...)
fracoes_vapor = np.linspace(0.0, 0.5, 11)

# Listas para guardar os resultados do final do reator (z = L)
H2_final = []
CO_final = []
razao_H2_CO = []

print("Iniciando a varredura da fração de vapor...")

for f_vap in fracoes_vapor:
    print(f"-> Testando fração de vapor: {f_vap*100:.1f}%")
    
    # Ajustar as frações molares da mistura de entrada
    f_ar = 1.0 - f_vap
    f_O2 = 0.21 * f_ar
    f_N2 = 0.79 * f_ar
    
    # [H2O, O2, CO, CO2, H2, CH4, N2, Tar]
    X_g0 = [f_vap, f_O2, 0.0, 0.0, 0.0, 0.0, f_N2, 0.0]
    
    # Recalcular as concentrações molares iniciais (C_g0)
    Mm_g0 = (materiais["Gases"]["H2O"]["MM"] * f_vap + 
             materiais["Gases"]["O2"]["MM"] * f_O2 + 
             materiais["Gases"]["N2"]["MM"] * f_N2) / 1000.0
    rhom_g0_real = (P0 * Mm_g0) / (Rg * T0)
    Cm_g0 = rhom_g0_real / Mm_g0
    C_g0 = [Cm_g0 * x for x in X_g0]
    
    # Montar vetor de estado inicial para a rodada atual
    y0 = [rho_B0, rho_char0, rho_M0] + C_g0 + [U_g0, P0, U_s0, T0, T0]
    
    # Resolver a EDO
    solucao = solve_ivp(
        fun=gasifier_odes,
        t_span=(0, reator_params["L"]),
        y0=y0,
        args=(materiais, reator_params),
        method='RK23',
        max_step=0.05
    )
    
    # Extrair os resultados da ÚLTIMA posição do reator (solucao.y[:,-1])
    # Os índices dos gases são: H2O (3), O2 (4), CO (5), CO2 (6), H2 (7)
    conc_CO_final = solucao.y[5][-1]
    conc_H2_final = solucao.y[7][-1]
    
    # Guardar nas listas
    CO_final.append(conc_CO_final)
    H2_final.append(conc_H2_final)
    
    # Calcular a razão H2 / CO
    razao = conc_H2_final / conc_CO_final if conc_CO_final > 0 else 0
    razao_H2_CO.append(razao)

print("Varredura concluída!")

# ==========================================
# 3. PLOTAGEM DOS RESULTADOS
# ==========================================
porcentagens_vapor = fracoes_vapor * 100

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Gráfico 1: Concentrações Finais
ax1.plot(porcentagens_vapor, H2_final, 'cyan', marker='o', label='H₂ Final', linewidth=2)
ax1.plot(porcentagens_vapor, CO_final, 'red', marker='s', label='CO Final', linewidth=2)
ax1.set_xlabel('Fração de Vapor Injetada (%)', fontsize=12)
ax1.set_ylabel('Concentração Final no Gaseificador (mol/m³)', fontsize=12)
ax1.set_title('Efeito da Injeção de Vapor na Produção de Gás', fontsize=14)
ax1.legend()
ax1.grid(True)

# Gráfico 2: Razão H2 / CO para Síntese de Metanol
ax2.plot(porcentagens_vapor, razao_H2_CO, 'purple', marker='^', linewidth=2)
ax2.axhline(y=2.0, color='green', linestyle='--', label='Razão Ideal para Metanol (2.0)')
ax2.set_xlabel('Fração de Vapor Injetada (%)', fontsize=12)
ax2.set_ylabel('Razão H₂ / CO', fontsize=12)
ax2.set_title('Adequação do Gás para o Reator de Metanol', fontsize=14)
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()