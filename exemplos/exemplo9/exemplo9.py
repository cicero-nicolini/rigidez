"""
EXEMPLO 9 - Exemplos do professor Masuero para cálculo por flexão simples

Objetivo: determinar armaduras

"""
from src.no import No
from src.concreto import Concreto
from src.aco import Aco
from src.secao import Secao
from src.barra import Barra

print("#########################################################################################################")
print("Teste Prova 1 2025/1 Questão 2")
print("#########################################################################################################")

#flexão simples, armadura dupla e concreto grupo 2

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(70)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("viga", 19, 50, concreto.Eci, d = 43, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
viga.fl[2] = 25000.0
viga.fl[5] = 25000.0

viga.fl[0] = 0.0
viga.fl[3] = 0.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Teste Exemplo 1 - aula 06b")
print("#########################################################################################################")

#flexão simples, armadura simples e concreto grupo 2

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(70)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("viga", 20, 40, concreto.Eci, d = 35, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
viga.fl[2] = 9000.0
viga.fl[5] = 9000.0

viga.fl[0] = 0.0
viga.fl[3] = 0.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)

print("#########################################################################################################")
print("Teste Exemplo 2 - aula 06b")
print("#########################################################################################################")

#flexão simples, armadura dupla e concreto grupo 2, unico parâmetro que muda é o esforço M

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
viga.fl[2] = 18000.0
viga.fl[5] = 18000.0

viga.fl[0] = 0.0
viga.fl[3] = 0.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(viga.ASL2)
print(viga.ASL)



