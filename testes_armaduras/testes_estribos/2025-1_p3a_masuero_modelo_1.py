import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rigidez import *

print("#########################################################################################################")
print("Teste Prova 3a Questão 2 - c)")
print("#########################################################################################################")

# estribos a 90, modelo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(35)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("viga", 19, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco, 1, 90.0)

# Esforços para teste
M = np.array([0.0, 0.0])
N = np.array([0.0, 0.0])
V = np.array([157.58, -108.67])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M, N)

# Passando esforços e parâmetros para metodo do objeto barra que calcula estribos e armazena no vetor ASW_S
viga.calcula_vetor_ASW_S(viga.modelo, viga.α, viga.θ, M, N, V)

# Visualização de resultados da area de estribo do objeto barra
print(viga.ASW_S)
