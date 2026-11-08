import matplotlib.pyplot as plt
from src.furo import Furo
from src.modelo_furos import ModeloFuros
from src.visualizacao import Visualizacao

furos = [
    Furo(lf=200, b=100, a=30, c=25),
    Furo(lf=500, b=120, a=25, c=30),
    Furo(lf=800, b=80, a=35, c=20),
]

cargas = [
    {"nome": "Peso próprio", "w": 0.15, "gamaf": 1.40},
    {"nome": "Revestimento", "w": 0.15, "gamaf": 1.40},
    {"nome": "Sobrecarga", "w": 0.15, "gamaf": 1.50},
]

modelo = ModeloFuros(
    bw=30.0,
    h=80.0,
    l=1100.0,
    furos=furos,
    ap1=40.0,
    ap2=40.0,
    cargas=cargas,
    fck=35.0,
    fyk=500.0,
    fywk=500.0,
    modelo=1,
    cob=3.5,
    espacamento_maximo=20.0,
    vinculacao_esquerda=1,
    vinculacao_direita=3
)

visualizacao = Visualizacao(modelo)

visualizacao.desenhar_geometria(
    mostrar_geometria=True,
    mostrar_modelo=False,
    mostrar_pontos_nos=False,
    mostrar_numeros_nos=False,
    mostrar_numeros_barras=False,
    mostrar_apoios=True,
    mostrar_cargas_distribuidas=False,
    mostrar_forcas_nodais=False,
    mostrar_momentos_nodais=False,
    mostrar_deformada=False,
    mostrar_momento=True,
    mostrar_normal=False,
    mostrar_cortante=False,
    mostrar_diagrama_as_longitudinal=False,
    mostrar_diagrama_as_transversal=False,
    mostrar_armadura_longitudinal=True,
    mostrar_armadura_transversal=True,
    mostrar_valores=True,
    mostrar_escalas=True,
    controles=True
)

plt.show()