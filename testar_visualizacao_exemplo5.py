"""Na raiz de rigidez: python testar_visualizacao_exemplo5.py"""
import matplotlib.pyplot as plt
from src.modelo import Modelo
from src.visualizacao import Visualizacao

# Mesmos parâmetros de exemplos/exemplo5/exemplo5.py
modelo = Modelo(
    bw=30.0, h=80.0, l=940.0, lf=150.0 + 60.0,
    a=25.0, b=120.0, c=30.0, ap1=40.0, ap2=40.0,
    w=0.45 + 0.06 + 0.02 - 0.0023,
    fck=35.0, fyk=500.0, fywk=500.0,
    modelo=1, α=90.0, θ=30.0, αE=1,espacamento_maximo = 20.0
)

visualizacao = Visualizacao(modelo)

visualizacao.desenhar_geometria(
    mostrar_nos=False,
    mostrar_barras=False,
    controles=True
)

plt.show()
