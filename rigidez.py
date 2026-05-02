import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

meu_nos = [
    (20,40),
    (40,40),
    (20,20),
    (40,20),
]

meu_barras = [
    (0,1),
    (1,3),
    (2,3),
    (0,2),
]

def desenha_estrutura(nos, barras):
    cord_x = []
    cord_y = []

    for no in nos:
        x, y = no
        cord_x.append(x)
        cord_y.append(y)

    linhas = []

    for conec in barras:
        indice_inicio, indice_fim = conec
        coord_inicio = nos[indice_inicio]
        coord_fim = nos[indice_fim]
        linhas.append([coord_inicio, coord_fim])

    lc = LineCollection(linhas, linewidths=2)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.plot(cord_x, cord_y, 'ro')
    plt.show()



desenha_estrutura(meu_nos, meu_barras)



