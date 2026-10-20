"""
EXEMPLO 4 - Exemplo referente ao primeiro teste da classe modelo

Objetivo: determinar deslocamentos, solicitações, reações e armaduras

"""
from src.modelo import Modelo

print("#########################################################################################################")
print("Teste Modelo")
print("#########################################################################################################")

modelo = Modelo(20.0, 50.0, 220.0, 100.0, 30.0, 40.0, 10.0, 20.0, 20.0, 10.0, 25.0, 500.0, 500.0, 1)
modelo.desenha_area_armaduras()
