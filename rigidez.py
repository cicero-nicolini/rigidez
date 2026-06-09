import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.collections import LineCollection


class No:
    def __init__(self, x, y, apoio = False):
        self.x = x
        self.y = y
        self.apoio = apoio

class Barra: 
    def __init__(self, no1, no2, cargas=[]):
        self.indice_no1 = no1
        self.indice_no2 = no2
        self.cargas = cargas

class CargaDistribuida: 
    def __init__(self, a, lw, w1, w2):
        '''
        a: distancia ate carga
        lw: comprimento carga
        w1: primeiro valor y da carga
        w2: segundo valor y
        '''
        self.a = a
        self.lw = lw
        self.w1 = w1
        self.w2 = w2

class CargaPontual: 
    def __init__(self, a, Px, Py):
        self.a = a
        self.Px = Px
        self.Py = Py

meu_nos = [
    No(60,135,True),
    No(160,135),
    No(260,60,True),
    #No(20,20),
]

c1 = CargaDistribuida(0, 100, 0.24, 0.24)
c2 = CargaPontual(62.5, 12, 16)

meu_barras = [
    Barra(0,1,[c1]),
    Barra(1,2,[c2]),
    #Barra(2,3),
    #Barra(0,2),
]

def comprimento_barra(no1, no2): #argamentos são coordenadas dos nós
    dx = no2.x - no1.x
    dy = no2.y - no1.y

    return math.sqrt( dx * dx + dy * dy )

def matriz_barra_local(E, A, I, L):
    kl = np.zeros((6,6))

    a1 = E*A/L
    a2 = (12.0*E*I)/L**3
    a3 = (6.0*E*I)/L**2
    a4 = (4.0*E*I)/L
    a5 = (2.0*E*I)/L

    kl[0,0] = a1                           
    kl[0,3] = -a1
    kl[1,1] = a2
    kl[1,2] = a3
    kl[1,4] = -a2
    kl[1,5] = a3
    kl[2,2] = a4
    kl[2,4] = -a3
    kl[2,5] = a5
    kl[3,3] = a1
    kl[4,4] = a2
    kl[4,5] = -a3
    kl[5,5] = a4

    kl[3,0] = kl[0,4]
    kl[2,1] = kl[1,2]
    kl[4,1] = kl[1,4]
    kl[5,1] = kl[1,5]
    kl[4,2] = kl[2,4]
    kl[5,2] = kl[2,5]
    kl[5,4] = kl[4,5]

    return kl

def matriz_rotacao(no1, no2): #argamentos são coordenadas dos nós
    r = np.zeros((6,6))

    dx = no2.x - no1.x
    dy = no2.y - no1.y
    l = math.sqrt( dx * dx + dy * dy )
    s = dy/l
    c = dx/l

    r[0,0] = c
    r[1,0] = -s
    r[0,1] = s
    r[1,1] = c
    r[2,2] = 1
    r[3,3] = c
    r[4,3] = -s
    r[3,4] = s
    r[4,4] = c
    r[5,5] = 1

    return r

def calcula_fepl_distribuido(a, lw, w1, w2, L):
    fepl = np.zeros((6))

    b = L - lw - a
    wm = (w1 + w2)/2.0
    wd = w2 - w1

    s1 = 10.0*((L*L + a*a)*(L+a)-(a*a+b*b)*(a-b)-(L*b)*(L+b)-a**3.0)
    s2 = lw*(L*(2.0*L+a+b)-3.0*(a-b)*(a-b)-2.0*a*b)
    s3 = 120.0*a*b*(a + lw)+10.0*lw*(6.0*a*a+4.0*L*lw-3.0*lw*lw)
    s4 = 10.0*L*lw*lw-10.0*lw*a*(L-3.0*b)-9.0*lw**3

    rb = (lw*(s1 * wm + s2 * wd ))/(20.0*L**3)
    ra = lw*wm-rb
    mb = -(lw*(s3*wm+s4*wd))/(120.0*L*L)
    ma = -mb+rb*L-a*lw*wm-(lw*lw*(2.0*w2+w1))/6.0

    fepl[1] = ra
    fepl[2] = ma
    fepl[4] = rb
    fepl[5] = mb

    return fepl

def calcula_fepl_pontual(a, Px, Py, L):
    fepl = np.zeros((6))

    b = L - a

    sa = L + 2*a
    sb = L + 2*b

    ha = -(Px*b*b*sa)/(L**3)
    hb = -(Px*a*a*sb)/(L**3)
    ra = (Py*b*b*sa)/(L**3)
    rb = (Py*a*a*sb)/(L**3)
    ma = (Py*a*b*b)/(L*L)
    mb = -(Py*a*a*b)/(L*L)

    fepl[0] = ha
    fepl[1] = ra
    fepl[2] = ma
    fepl[3] = hb
    fepl[4] = rb
    fepl[5] = mb

    return fepl

def calcula_q(indice1, indice2): #correspondencia graus de liberdade da barra para estrutura
    q = [0]*6
    z = 0
    M = np.array([indice1, indice2])
    for j in range(0,2):
        for jk in range(0,3):
            q[z] = 3*(M[j])+jk
            z += 1

    return q

def remove_grau_liberdade(k_estrutura, grau_liberdade):
    tamanho = k_estrutura.shape[0]
    for i in range(tamanho):
        k_estrutura[grau_liberdade, i] = 0.0
        k_estrutura[i, grau_liberdade] = 0.0
    k_estrutura[grau_liberdade,grau_liberdade] = 1.0


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
dimensao = len(meu_nos)*3
k_estrutura = np.zeros((dimensao, dimensao))
forcas = np.zeros(dimensao)
graus_restritos = []

for barra in meu_barras:
    no_inicio = meu_nos[barra.indice_no1] #coordenadas primeiro nó da barra
    no_fim = meu_nos[barra.indice_no2] #coordenadas segundo nó da barra
    R = matriz_rotacao(no_inicio, no_fim)
    R_T = R.T  
    E = 10000.0
    area = 2.0 * 5.0
    inercia = 1000.0
    tamanho = comprimento_barra(no_inicio, no_fim)
    KL = matriz_barra_local(E, area, inercia, tamanho)
    K = R_T @ KL @ R

    vetorq = calcula_q(barra.indice_no1, barra.indice_no2)
    
    if no_inicio.apoio:
        graus_restritos.append(vetorq[0])
        graus_restritos.append(vetorq[1])
        graus_restritos.append(vetorq[2])
    if no_fim.apoio:
        graus_restritos.append(vetorq[3])
        graus_restritos.append(vetorq[4])
        graus_restritos.append(vetorq[5])
    
    fep = np.zeros((6))
    for carga in barra.cargas:
        if isinstance(carga, CargaDistribuida):
            fep += calcula_fepl_distribuido(carga.a, carga.lw, carga.w1, carga.w2, tamanho)
        else:
            fep += calcula_fepl_pontual(carga.a, carga.Px, carga.Py, tamanho)
    fep = R_T @ fep
    print(fep)
    for i in range(len(fep)):
        forcas[vetorq[i]] = fep[i]

    for j in range(0,6):
        for jk in range(0,6):
            k_estrutura[vetorq[j], vetorq[jk]] += K[j,jk]


for grau in graus_restritos:
    remove_grau_liberdade(k_estrutura, grau)
    forcas[grau] = 0.0


print(k_estrutura)
print(forcas)

deslocamentos = np.linalg.solve(k_estrutura, forcas)
print(deslocamentos)
  







    





    




