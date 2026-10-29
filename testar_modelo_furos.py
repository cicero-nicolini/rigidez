import matplotlib.pyplot as plt
from src.furo import Furo
from src.modelo_furos import ModeloFuros
from src.visualizacao import Visualizacao

furos = [
    Furo(lf=200, b=100, a=30, c=25),
    Furo(lf=500, b=120, a=25, c=30),
    Furo(lf=800, b=80, a=35, c=20),
]

modelo = ModeloFuros(
    bw=30.0,
    h=80.0,
    l=1100.0,
    furos=furos,
    ap1=40.0,
    ap2=40.0,
    w=0.45,
    fck=35.0,
    fyk=500.0,
    fywk=500.0,
    modelo=1,
    espacamento_maximo=20.0
)

visualizacao = Visualizacao(modelo)

visualizacao.desenhar_geometria(
    mostrar_nos=False,
    mostrar_barras=False,
    controles=True
)

plt.show()