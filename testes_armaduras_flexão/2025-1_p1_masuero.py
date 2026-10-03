from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Prova 1 Questão 2")
print("#########################################################################################################")

#flexão simples, armadura dupla e concreto grupo 2

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(70)
aco = Aco(500)

# Definição da seção transversal
secao = Secao("viga", 19, 50, concreto.Eci, d = 43, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M1 = np.array([25000.0, 0.0])
N1 = np.array([0, 0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M1, N1)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)