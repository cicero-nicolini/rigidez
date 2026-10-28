"""
EXEMPLO 8 - Exemplos do professor Masuero para cálculo de estribos modelo 2

Objetivo: determinar armaduras

"""
from src.no import No
from src.concreto import Concreto
from src.aco import Aco
from src.secao import Secao
from src.barra import Barra

print("#########################################################################################################")
print("Exemplo 1 - b) aula sobre dimensionamento ao cortante - modelo 2")
print("#########################################################################################################")

# estribos a 90, modelo 2 - Masuero

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(30)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("viga", 18, 45, concreto.Eci, d = 40, d_linha = 5)

# Definição da barra para obtencao de resultados
viga = Barra(1, no1, no2, secao, concreto, aco, 2, 90.0, 45.0)

# Esforços para teste
viga.fl[1] = 105.0
viga.fl[4] = 105.0

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




