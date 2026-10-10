from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Exemplo 1")
print("#########################################################################################################")

#flexo-tração com grande excentricidade, armadura simples e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500)

# Definição da seção transversal
secao = Secao("viga", 20, 60, concreto.Eci, d = 55, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M = np.array([20000.0])
N = np.array([500.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M, N)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Exemplo 2")
print("#########################################################################################################")

#O teste bate com o exemplo apenas sem considerar
#xlim de acordo com o limite de ductilidade da NBR6118:2026

#flexo-tração com grande excentricidade, armadura dupla e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500)

# Definição da seção transversal
secao = Secao("viga", 20, 60, concreto.Eci, d = 55, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M = np.array([50000.0])
N = np.array([600.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M, N)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Exemplo 3")
print("#########################################################################################################")

#flexo-tração com pequena excentricidade e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500)

# Definição da seção transversal
secao = Secao("viga", 20, 60, concreto.Eci, d = 55, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M = np.array([10000.0])
N = np.array([500.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M, N)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)