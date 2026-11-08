"""Na raiz de rigidez: python testar_visualizacao_exemplo5.py"""
import matplotlib.pyplot as plt
from src.modelo import Modelo
from src.visualizacao import Visualizacao

# Mesmos parâmetros de exemplos/exemplo5/exemplo5.py
modelo = Modelo(
    bw=30.0, 
    h=80.0, 
    l=940.0, 
    lf=150.0 + 60.0,
    a=25.0, 
    b=120.0, 
    c=30.0, 
    ap1=30.0, 
    ap2=40.0,
    w=0.45 + 0.06 + 0.02 - 0.0023,
    fck=35.0, 
    fyk=500.0, 
    fywk=500.0,
    modelo=1, 
    cob=3.5,
    α=90.0, 
    θ=30.0, 
    αE=1,espacamento_maximo = 20.0,
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
    mostrar_diagrama_decalado=True,
    decalagem_automatica=True,
    ancoragem_automatica=True,
    boa_aderencia=True,
    ancoragem_apoios=True,
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
