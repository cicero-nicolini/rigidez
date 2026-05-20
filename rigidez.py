import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.collections import LineCollection


class No:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Barra: 
    def __init__(self, no1, no2):
        self.indice_no1 = no1
        self.indice_no2 = no2        

meu_nos = [
    No(10,40),
    No(30,40),
    No(15,20),
    No(20,20),
]

meu_barras = [
    Barra(0,1),
    Barra(1,3),
    Barra(2,3),
    Barra(0,2),
]

def desenha_estrutura(nos, barras):
    cord_x = [] 
    cord_y = [] 

    for no in nos:
        cord_x.append(no.x) 
        cord_y.append(no.y) 

    linhas = []

    for barra in barras:
        no_inicio = nos[barra.indice_no1] #coordenadas primeiro nó da barra
        no_fim = nos[barra.indice_no2] #coordenadas segundo nó da barra
        no1 = (no_inicio.x, no_inicio.y)
        no2 = (no_fim.x, no_fim.y)
        linhas.append([no1, no2])

    lc = LineCollection(linhas, linewidths=2)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.plot(cord_x, cord_y, 'ro')
    plt.show()



desenha_estrutura(meu_nos, meu_barras)

def comprimento_barra(no1, no2):
    dx = no2.x - no1.x
    dy = no2.y - no1.y

    return math.sqrt( dx * dx + dy * dy )


for barra in meu_barras:
    no_inicio = meu_nos[barra.indice_no1] #coordenadas primeiro nó da barra
    no_fim = meu_nos[barra.indice_no2] #coordenadas segundo nó da barra 
    tamanho = comprimento_barra(no_inicio, no_fim)
    print(tamanho)




    





    




