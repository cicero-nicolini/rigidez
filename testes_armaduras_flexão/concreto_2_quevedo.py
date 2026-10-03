from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Prova 1 Exemplo 3")
print("#########################################################################################################")

#flexo-tração com grande excentricidade, armadura simples e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(36)
aco = Aco(250)

# Definição da seção transversal
secao = Secao("viga", 35, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
M1 = np.array([30000.0, 0.0])
N1 = np.array([170.0, 0.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL(M1, N1)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Teste Prova 2 Questão 1")
print("#########################################################################################################")

#O teste funciona apenas sem considerar
#xlim de acordo com o limite de ductilidade da NBR6118:2026

#flexo-compressão com grande excentricidade, armadura dupla e concreto grupo 1

# Definição dos materiais para uma seção de concreto armado
concreto2 = Concreto(40)
aco2 = Aco(400)

# Definição da seção transversal
secao2 = Secao("viga", 20, 80, concreto2.Eci, d = 75, d_linha = 5)

# Definição da barra para obtencao de resultados
viga2 = Barra(2, no1, no2, secao2, concreto2, aco2)

# Esforços para teste
M2 = np.array([80000.0, 0.0])
N2 = np.array([-2400.0, 0.0])

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga2.calcula_vetor_ASL(M2, N2)

# Visualização de resultados das areas de armaduras do objeto barra
print(viga2.ASL2)
print(viga2.ASL)
