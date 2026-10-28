"""
EXEMPLO 7 - Exemplos do professor Masuero para cálculo de estribos modelo 1

Objetivo: determinar armaduras

"""
from src.no import No
from src.concreto import Concreto
from src.aco import Aco
from src.secao import Secao
from src.barra import Barra

print("#########################################################################################################")
print("Teste Prova 3a 2025/1 Questão 2 - c)")
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
viga.fl[1] = 157.58
viga.fl[4] = -108.67

viga.fl[2] = 0.0
viga.fl[5] = 0.0

viga.fl[0] = 0.0
viga.fl[3] = 0.0

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
viga.calcula_vetor_ASL()

# Passando esforços e parâmetros para metodo do objeto barra que calcula estribos e armazena no vetor ASW_S
viga.calcula_vetor_ASW_S(viga.modelo, viga.α, viga.θ)

# Visualização de resultados da area de estribo do objeto barra
print(viga.ASW_S)







