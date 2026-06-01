import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.collections import LineCollection


class No:
    def __init__(self, num, x, y):
        self.num = num 
        self.x = x
        self.y = y

class Barra: 
    def __init__(self, num, no1, no2,E,A,I):
        self.num = num
        self.noi = no1
        self.noj = no2
        self.E = E
        self.A = A
        self.I = I
        self.L = None      
        self.kl = np.zeros((6,6))
        self.r = np.zeros((6,6))
        self.fepl = np.zeros((6))

    def comprimento_barra(self):
        dx = self.noj.x - self.noi.x
        dy = self.noj.y - self.noi.y
        self.L = math.sqrt( dx * dx + dy * dy )

    def calcula_klocal(self):
        
        a1 = self.E*self.A/self.L
        a2 = (12.0*self.E*self.I)/self.L**3
        a3 = (6.0*self.E*self.I)/self.L**2
        a4 = (4.0*self.E*self.I)/self.L
        a5 = (2.0*self.E*self.I)/self.L

        self.kl[0,0] = a1
        self.kl[0,4] = -a1
        self.kl[1,1] = a2
        self.kl[1,2] = a3
        self.kl[1,4] = -a2
        self.kl[1,5] = a3
        self.kl[2,2] = a4
        self.kl[2,4] = -a3
        self.kl[2,5] = a5
        self.kl[3,3] = a1
        self.kl[4,4] = a2
        self.kl[4,5] = -a3
        self.kl[5,5] = a4

        self.kl[4,0] = self.kl[0,4]
        self.kl[2,1] = self.kl[1,2]
        self.kl[4,1] = self.kl[1,4]
        self.kl[5,1] = self.kl[1,5]
        self.kl[4,2] = self.kl[2,4]
        self.kl[5,2] = self.kl[2,5]
        self.kl[5,4] = self.kl[4,5]        

    def calcula_r(self): 
        dx = self.noj.x - self.noi.x
        dy = self.noj.y - self.noi.y
        s = dy/self.L
        c = dx/self.L

        self.r[0,0] = c
        self.r[1,0] = -s
        self.r[0,1] = s
        self.r[1,1] = c   
        self.r[2,2] = 1
        self.r[3,3] = c
        self.r[4,3] = -s
        self.r[3,4] = s      
        self.r[4,4] = c
        self.r[5,5] = 1 
    

class Carregamento_distribuido: 
     
    def __init__(self, a, lw, w1, w2, barra):
        self.a = a
        self.lw = lw
        self.w1 = w1
        self.w2 = w2
        self.barra = barra
        
    def calcula_Fepl(self):
        L = self.barra.L
        b = L - self.lw - self.a
        wm = (self.w1 + self.w2)/2.0
        wd = self.w2 - self.w1
        s1 = 10.0*((L*L + self.a*self.a)*(L+self.a)-(self.a*self.a+b*b)*(self.a-b)-(L*b)*(L+b)-self.a**3.0)
        s2 = self.lw*(L*(2.0*L+self.a+b)-3.0*(self.a-b)*(self.a-b)-2.0*self.a*b)
        s3 = 120.0*self.a*b*(self.a + self.lw)+10.0*self.lw*(6.0*self.a*self.a+4.0*L*self.lw-3.0*self.lw*self.lw)
        s4 = 10.0*L*self.lw*self.lw-10.0*self.lw*self.a*(L-3.0*b)-9.0*self.lw**3
        rb = (self.lw*(s1 * wm + s2 * wd ))/(20.0*L**3)
        ra = self.lw*wm-rb
        mb = -(self.lw*(s3*wm+s4*wd))/(120.0*L*L)
        ma = -mb+rb*L-self.a*self.lw*wm-(self.lw*self.lw*(2.0*self.w2+self.w1))/6.0
        self.barra.fepl[1] = ra
        self.barra.fepl[2] = ma
        self.barra.fepl[4] = rb
        self.barra.fepl[5] = mb
        #self.barra.fepl[0] = self.barra.fepl[0]+ 1 

class Carregamento_pontual: 
     
    def __init__(self, a, P, barra):
        self.a = a
        self.P = P
        self.barra = barra

    def calcula_Fepl(self):
        L = self.barra.L
        b = L - self.a
        sa = L + 2*self.a
        sb = L + 2*b
        ra = (self.P*b*b*sa)/(L**3)
        rb = (self.P*self.a*self.a*sb)/(L**3)
        ma = (self.P*self.a*b*b)/(L*L)
        mb = -(self.P*self.a*self.a*b)/(L*L)
        self.barra.fepl[1] = ra
        self.barra.fepl[2] = ma
        self.barra.fepl[4] = rb
        self.barra.fepl[5] = mb    
        

no1 = No(1,60.0,135.0)
no2 = No(2,160.0,135.0)
no3 = No(3,260.0,60.0)

barra1 = Barra(1.0,no1,no2,10000.0,2.0*5.0,1000.0)
barra2 = Barra(2.0,no2,no3,10000.0,2.0*5.0,1000.0)
barra1.comprimento_barra()
barra2.comprimento_barra()
barra1.calcula_klocal()
barra1.calcula_klocal()
barra1.calcula_r()
print(barra1.kl)
print(barra1.r)
carregamento1 = Carregamento_distribuido(0,100,0.24,0.24,barra1)
carregamento2 = Carregamento_pontual(62.5,16.0,barra2)
carregamento1.calcula_Fepl()
carregamento2.calcula_Fepl()
print(barra1.fepl)
print(barra2.fepl)
#carregamento1.mostre()
#carregamento1.mostre()
#carregamento1.mostre()

#    E = 10000.0
#    area = 2.0 * 5.0
#    inercia = 1000.0


#class Estrutura:


#meu_nos = [
#    No(60,135),
#    No(160,135),
#    No(260,60),
#    #No(20,20),
#]
#
#meu_barras = [
#    Barra(0,1),
#    Barra(1,2),
#    #Barra(2,3),
#    #Barra(0,2),
#]
#
##cargas = []
#
#def comprimento_barra(no1, no2): #argamentos são coordenadas dos nós
#    dx = no2.x - no1.x
#    dy = no2.y - no1.y
#
#    return math.sqrt( dx * dx + dy * dy )
#
#def momento_inercia(b, l):
#
#    return b * l**3/12
#
#def matriz_local(E, A, I, L):
#
#    k = np.array([
#        [ (E*A)/L,            0.0,                 0.0,  -(E*A)/L,               0.0,                 0.0],
#
#        [     0.0,      (12.0*E*I)/L**3,  6.0*E*I/L**2,       0.0,  -(12.0*E*I)/L**3,      (6.0*E*I)/L**2],
#
#        [     0.0,       (6.0*E*I)/L**2,   (4.0*E*I)/L,       0.0,   -(6.0*E*I)/L**2,         (2.0*E*I)/L],
#
#        [-(E*A)/L,                  0.0,           0.0,   (E*A)/L,               0.0,                 0.0],
#
#        [     0.0,     -(12.0*E*I)/L**3,   -6*E*I/L**2,       0.0,   (12.0*E*I)/L**3,     -(6.0*E*I)/L**2],
#
#        [     0.0,       (6.0*E*I)/L**2,       2*E*I/L,       0.0,   -(6.0*E*I)/L**2,         (4.0*E*I)/L]
#    ])
#
#    return k
#
#def matriz_R(no1, no2): #argamentos são coordenadas dos nós
#
#    dx = no2.x - no1.x
#    dy = no2.y - no1.y
#    l = math.sqrt( dx * dx + dy * dy )
#    seno = dy/l
#    cos = dx/l
#
#    r = np.array([
#        [  cos,             seno,           0.0,             0.0,            0.0,           0.0],
#
#        [-seno,              cos,           0.0,             0.0,            0.0,           0.0],
#
#        [  0.0,              0.0,           1.0,             0.0,            0.0,           0.0],
#
#        [  0.0,              0.0,           0.0,             cos,           seno,           0.0],
#
#        [  0.0,              0.0,           0.0,           -seno,            cos,           0.0],
#
#        [  0.0,              0.0,           0.0,             0.0,            0.0,           1.0]
#    ])
#
#    return r
#
#def Fep_local(tipo, valor, no1, no2): 
#
#    """
#    argumentos:
#
#    tipo: distribuida ou pontual
#    coordenadas dos nós
#
#    """
#    dx = no2.x - no1.x
#    dy = no2.y - no1.y
#    l = math.sqrt( dx * dx + dy * dy )
#    seno = dy/l
#    cos = dx/l
#
#    if tipo == "distribuida":
#        qx = valor*seno
#        qy = valor*cos 
#        ha = (qx*l)/2.0
#        va = (qy*l)/2.0
#        ma = (qy*l*l)/12.0
#        hb = (qx*l)/2.0
#        vb = (qy*l)/2.0
#        mb = -(qy*l*l)/12.0    
#        Fep_local = np.array([ha, va, ma, hb, vb, mb])
#        return Fep_local
#
#    elif tipo == "pontual":
#        px = valor*seno
#        py = valor*cos
#        ha = -px/2.0
#        va =  py/2.0
#        ma =  (py*l)/8.0
#        hb = -px/2.0
#        vb =  py/2.0
#        mb = -(py*l)/8.0    
#        Fep_local = np.array([ha, va, ma, hb, vb, mb])
#        return Fep_local
#
#def desenha_estrutura(nos, barras):
#    cord_x = [] 
#    cord_y = [] 
#
#    for no in nos:
#        cord_x.append(no.x) 
#        cord_y.append(no.y) 
#
#    linhas = []
#
#    for barra in barras:
#        no_inicio = nos[barra.indice_no1] #coordenadas primeiro nó da barra
#        no_fim = nos[barra.indice_no2] #coordenadas segundo nó da barra
#        no1 = (no_inicio.x, no_inicio.y)
#        no2 = (no_fim.x, no_fim.y)
#        linhas.append([no1, no2])
#
#    lc = LineCollection(linhas, linewidths=2)
#    fig, ax = plt.subplots()
#    ax.add_collection(lc)
#    ax.plot(cord_x, cord_y, 'ro')
#    plt.show()
#
#desenha_estrutura(meu_nos, meu_barras)
#
#for barra in meu_barras:
#    no_inicio = meu_nos[barra.indice_no1] #coordenadas primeiro nó da barra
#    no_fim = meu_nos[barra.indice_no2] #coordenadas segundo nó da barra
#    R = matriz_R(no_inicio, no_fim)
#    R_T = np.linalg.inv(R)   
#    E = 10000.0
#    area = 2.0 * 5.0
#    inercia = 1000.0
#    tamanho = comprimento_barra(no_inicio, no_fim)
#    KL = matriz_local(E, area, inercia, tamanho)
#    K = R_T @ KL @ R
#    print(tamanho)
#
#tipo = "distribuida"
#fep = Fep_local(tipo, 0.24, no_inicio, no_fim)
#print(fep)



    





    




