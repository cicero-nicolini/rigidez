from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Exemplo 1")
print("#########################################################################################################")

#flexão simples, armadura simples e concreto grupo 2

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(70)
aco = Aco(500)

# Definição da seção transversal
secao = Secao("viga", 20, 40, concreto.Eci, d = 35, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M1 = np.array([9000.0, 0.0])
N1 = np.array([0, 0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M1, N1)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Teste Exemplo 2")
print("#########################################################################################################")

#flexão simples, armadura dupla e concreto grupo 2, unico parâmetro que muda é o esforço M

# Definição da barra para obtencao de resultados
viga2 = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M2 = np.array([18000.0, 0.0])
N2 = np.array([0.0, 0.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga2.calcula_vetor_ASL(M2, N2)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga2.ASL2)
print(viga2.ASL)