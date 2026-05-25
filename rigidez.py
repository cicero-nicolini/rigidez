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

#class Estrutura:
 


meu_nos = [
    No(60,135),
    No(160,135),
    No(260,60),
    #No(20,20),
]

meu_barras = [
    Barra(0,1),
    Barra(1,2),
    #Barra(2,3),
    #Barra(0,2),
]

#cargas = []

def comprimento_barra(no1, no2): #argamentos são coordenadas dos nós
    dx = no2.x - no1.x
    dy = no2.y - no1.y

    return math.sqrt( dx * dx + dy * dy )

def momento_inercia(b, l):

    return b * l**3/12

def matriz_local(E, A, I, L):

    k = np.array([
        [ (E*A)/L,            0.0,                 0.0,  -(E*A)/L,               0.0,                 0.0],

        [     0.0,      (12.0*E*I)/L**3,  6.0*E*I/L**2,       0.0,  -(12.0*E*I)/L**3,      (6.0*E*I)/L**2],

        [     0.0,       (6.0*E*I)/L**2,   (4.0*E*I)/L,       0.0,   -(6.0*E*I)/L**2,         (2.0*E*I)/L],

        [-(E*A)/L,                  0.0,           0.0,   (E*A)/L,               0.0,                 0.0],

        [     0.0,     -(12.0*E*I)/L**3,   -6*E*I/L**2,       0.0,   (12.0*E*I)/L**3,     -(6.0*E*I)/L**2],

        [     0.0,       (6.0*E*I)/L**2,       2*E*I/L,       0.0,   -(6.0*E*I)/L**2,         (4.0*E*I)/L]
    ])

    return k

def matriz_R(no1, no2): #argamentos são coordenadas dos nós

    dx = no2.x - no1.x
    dy = no2.y - no1.y
    l = math.sqrt( dx * dx + dy * dy )
    seno = dy/l
    cos = dx/l

    r = np.array([
        [  cos,             seno,           0.0,             0.0,            0.0,           0.0],

        [-seno,              cos,           0.0,             0.0,            0.0,           0.0],

        [  0.0,              0.0,           1.0,             0.0,            0.0,           0.0],

        [  0.0,              0.0,           0.0,             cos,           seno,           0.0],

        [  0.0,              0.0,           0.0,           -seno,            cos,           0.0],

        [  0.0,              0.0,           0.0,             0.0,            0.0,           1.0]
    ])

    return r

def Fep_local(tipo, valor, no1, no2): 

    """
    argumentos:

    tipo: distribuida ou pontual
    coordenadas dos nós

    """
    dx = no2.x - no1.x
    dy = no2.y - no1.y
    l = math.sqrt( dx * dx + dy * dy )
    seno = dy/l
    cos = dx/l

    if tipo == "distribuida":
        qx = valor*seno
        qy = valor*cos 
        ha = (qx*l)/2.0
        va = (qy*l)/2.0
        ma = (qy*l*l)/12.0
        hb = (qx*l)/2.0
        vb = (qy*l)/2.0
        mb = -(qy*l*l)/12.0    
        Fep_local = np.array([ha, va, ma, hb, vb, mb])
        return Fep_local

    elif tipo == "pontual":
        px = valor*seno
        py = valor*cos
        ha = -px/2.0
        va =  py/2.0
        ma =  (py*l)/8.0
        hb = -px/2.0
        vb =  py/2.0
        mb = -(py*l)/8.0    
        Fep_local = np.array([ha, va, ma, hb, vb, mb])
        return Fep_local

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

for barra in meu_barras:
    no_inicio = meu_nos[barra.indice_no1] #coordenadas primeiro nó da barra
    no_fim = meu_nos[barra.indice_no2] #coordenadas segundo nó da barra
    R = matriz_R(no_inicio, no_fim)
    R_T = np.linalg.inv(R)   
    E = 10000.0
    area = 2.0 * 5.0
    inercia = 1000.0
    tamanho = comprimento_barra(no_inicio, no_fim)
    KL = matriz_local(E, area, inercia, tamanho)
    K = R_T @ KL @ R
    print(tamanho)

tipo = "distribuida"
fep = Fep_local(tipo, 0.24, no_inicio, no_fim)
print(fep)



    





    




