"""
EXEMPLO 10 - Exemplos do professor Américo para cálculo de armaduras por flexo-compressão

Objetivo: determinar armaduras

"""
from src.no import No
from src.concreto import Concreto
from src.aco import Aco
from src.secao import Secao
from src.barra import Barra

print("#########################################################################################################")
print("Exemplo 1")
print("#########################################################################################################")

#flexo-compressão com grande excentricidade, armadura simples e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("pilar", 25, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
pilar = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
pilar.fl[2] = 7000.0
pilar.fl[5] = 7000.0

pilar.fl[0] = 100.0
pilar.fl[3] = -100.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
pilar.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(pilar.ASL2)
print(pilar.ASL)

print("#########################################################################################################")
print("Exemplo 2")
print("#########################################################################################################")

#O teste bate com o exemplo apenas sem considerar
#xlim de acordo com o limite de ductilidade da NBR6118:2026

#flexo-compressão com grande excentricidade, armadura dupla e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("pilar", 25, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
pilar = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
pilar.fl[2] = 15000.0
pilar.fl[5] = 15000.0

pilar.fl[0] = 800.0
pilar.fl[3] = -800.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
pilar.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(pilar.ASL2)
print(pilar.ASL)

print("#########################################################################################################")
print("Exemplo 3")
print("#########################################################################################################")

#flexo-compressão, não há necessidade teórica de armadura e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("pilar", 25, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
pilar = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
pilar.fl[2] = 3000.0
pilar.fl[5] = 3000.0

pilar.fl[0] = 630.0
pilar.fl[3] = -630.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
pilar.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(pilar.ASL2)
print(pilar.ASL)

print("#########################################################################################################")
print("Exemplo 4")
print("#########################################################################################################")

#flexo-compressão com pequena excentricidade e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("pilar", 25, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
pilar = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
pilar.fl[2] = 10000.0
pilar.fl[5] = 10000.0

pilar.fl[0] = 1250.0
pilar.fl[3] = -1250.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
pilar.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(pilar.ASL2)
print(pilar.ASL)

print("#########################################################################################################")
print("Exemplo 5")
print("#########################################################################################################")

#compressão composta e concreto grupo 1

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(25)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("pilar", 25, 50, concreto.Eci, d = 45, d_linha = 5)

# Definição da barra para obtencao de resultados
pilar = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
pilar.fl[2] = 10000.0
pilar.fl[5] = 10000.0

pilar.fl[0] = 2000.0
pilar.fl[3] = -2000.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
pilar.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(pilar.ASL2)
print(pilar.ASL)




