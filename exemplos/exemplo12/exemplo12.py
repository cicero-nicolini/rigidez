"""
EXEMPLO_SCH - Exemplo referente ao material do professor Schwingel

Objetivo: determinar armaduras

"""
from src.no import No
from src.concreto import Concreto
from src.aco import Aco
from src.secao import Secao
from src.barra import Barra

print("#########################################################################################################")
print("Modelo Schwingel Banzo Inferior")
print("#########################################################################################################")

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)

# Definição dos materiais para uma seção de concreto armado
concreto = Concreto(35)
aco = Aco(500, 500)

# Definição da seção transversal
secao = Secao("banzo inferior", 30, 30, concreto.Eci, d = 26.5, d_linha = 3.5)

# Definição da barra para obtencao de resultados
banzo = Barra(1, no1, no2, secao, concreto, aco)

# Esforços para teste
banzo.fl[2] = 6611.0
banzo.fl[5] = 6611.0

banzo.fl[0] = -735.36
banzo.fl[3] = 735.36

# Passando esforços para metodo do objeto barra que calcula armadura e armazena no vetor ASL e ASL2
banzo.calcula_vetor_ASL()

# Visualização de resultados das areas de armaduras do objeto barra
print(banzo.ASL2)
print(banzo.ASL)




