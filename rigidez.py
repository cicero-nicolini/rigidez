import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.collections import LineCollection


class No:
    def __init__(self, num, x, y):
        self.num = int(num) 
        self.x = x
        self.y = y
        self.Fx = 0.0
        self.Fy = 0.0
        self.Mz = 0.0
        self.Tx = False
        self.Ty = False
        self.Rz = False

class Barra: 
    def __init__(self, num, no1, no2,E,A,I):
        self.num = int(num)
        self.noi = no1
        self.noj = no2
        self.E = E
        self.A = A
        self.I = I
        self.L = None      
        self.kl = np.zeros((6,6))
        self.k = np.zeros((6,6))
        self.r = np.zeros((6,6))
        self.fepl = np.zeros((6))
        self.fep = np.zeros((6))
        self.q = np.zeros((6),dtype=int)

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
        self.kl[0,3] = -a1
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

        self.kl[3,0] = self.kl[0,3]
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

    def calcula_k(self):
        self.k = np.linalg.inv(self.r)@ self.kl @ self.r
    
    def calcula_q(self):
        z = -1
        M = np.array([self.noi.num, self.noj.num])
        for j in range(0,2):
            for jk in range(0,3):
                z = z + 1
                self.q[z] = 3*(M[j]-1)+jk
        
        
        
class Estrutura:
    def __init__(self,nos,barras):
        self.nos = nos
        self.barras = barras
        self.nnos = len(nos)
        self.fnos = np.zeros((self.nnos*3))
        self.f01 = np.zeros((self.nnos*3))
        self.k = np.zeros((self.nnos*3,self.nnos*3), dtype=int)
        self.k01 = np.zeros((self.nnos*3,self.nnos*3), dtype=int)
                
    def monta_k(self):
        for barra in self.barras:
            barra.comprimento_barra()
            barra.calcula_klocal()
            barra.calcula_r()
            barra.calcula_k()
            barra.calcula_q()
            for j in range(0,6):
                for jk in range(0,6):
                    self.k[barra.q[j], barra.q[jk]] = self.k[barra.q[j], barra.q[jk]]+barra.k[j,jk]

        #print(self.k)

    def monta_k01(self):
        self.k01 = self.k
        for no in self.nos:
            if no.Tx == True:
                gdl = 3*no.num-3
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0
            if no.Ty == True:
                gdl = 3*no.num-2
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0               
            if no.Rz == True:
                gdl = 3*no.num-1
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self. k01[gdl,gdl] = 1.0           
        
    def monta_fnos(self):    
        for no in self.nos:
            self.fnos[3*no.num-3] = no.Fx
            self.fnos[3*no.num-2] = no.Fy
            self.fnos[3*no.num-1] = no.Mz 
            
        for barra in self.barras:
             for j in range(0,6):
                self.fnos[barra.q[j]] += -barra.fep[j]
                
    def aplica_cc_fnos(self):
        for no in self.nos:
            if no.Tx == True:
                self.fnos[3*no.num-3] = 0
            if no.Ty == True:
                self.fnos[3*no.num-2] = 0
            if no.Rz == True:
                 self.fnos[3*no.num-1] = 0                                                
class Carregamento_distribuido: 
     
    def __init__(self, a, lw, w1, w2, barra):
        self.a = a
        self.lw = lw
        self.w1 = w1
        self.w2 = w2
        self.barra = barra
        
    def calcula_fepl(self):
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
        self.barra.fepl[1] += ra
        self.barra.fepl[2] += ma
        self.barra.fepl[4] += rb
        self.barra.fepl[5] += mb
        #self.barra.fepl[0] = self.barra.fepl[0]+ 1

    def calcula_fep(self):
        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl

class Carregamento_pontual:

    '''
        Py: componente y da força
        Px componente x da força

    ''' 
     
    def __init__(self, a, Px, Py, barra):
        self.a = a
        self.Py = Py
        self.Px = Px
        self.barra = barra

    def calcula_fepl(self):
        L = self.barra.L
        b = L - self.a
        sa = L + 2*self.a
        sb = L + 2*b

        ha = -(self.Px*b*b*sa)/(L**3)
        hb = -(self.Px*self.a*self.a*sb)/(L**3)
        ra = (self.Py*b*b*sa)/(L**3)
        rb = (self.Py*self.a*self.a*sb)/(L**3)
        ma = (self.Py*self.a*b*b)/(L*L)
        mb = -(self.Py*self.a*self.a*b)/(L*L)

        self.barra.fepl[0] += ha
        self.barra.fepl[1] += ra
        self.barra.fepl[2] += ma
        self.barra.fepl[3] += hb
        self.barra.fepl[4] += rb
        self.barra.fepl[5] += mb

    def calcula_fep(self):
        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl



# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)
no3 = No(3,200.0,0.0)
nos = [no1, no2, no3]

# Aplicação das cargas nodais
no2.Fy = -10
no2.Mz = -1000

# Aplicação das restrições nodais
no1.Tx = True
no1.Ty = True
no1.Rz = True
no3.Tx = True
no3.Ty = True
no3.Rz = True

# Definição das barras e propriedades
barra1 = Barra(1.0,no1,no2,10000.0,2.0*5.0,1000.0)
barra2 = Barra(2.0,no2,no3,10000.0,2.0*5.0,1000.0)
barras = [barra1, barra2]

# Definição dos carregamentos nas barras
carregamento1 = Carregamento_distribuido(0,100,0.24,0.24,barra1)
carregamento2 = Carregamento_pontual(62.5,12.0,16.0,barra2)

# Definição da estrutura
portico = Estrutura(nos,barras,grestritos)

# Monta matriz de rigidez global
portico.monta_k()
print(portico.k)

# Monta matriz de rigidez global com as condições de contorno
portico.monta_k01()
print(portico.k01)

# Calcula as forças de engastamento perfeito
carregamento1.calcula_fepl()
carregamento1.calcula_fep()
carregamento2.calcula_fepl()
carregamento2.calcula_fep()

# Monta vetor de cargas nodais global
portico.monta_fnos()

# Aplica condições de contorno no vetor de cargas nodais global
portico.aplica_cc_fnos()


print(portico.fnos)

#forcas_nos.calcula_fnos()
#forcas_nos.calcula_f01()

#print(barra1.fepl)
#print(barra2.fepl)
##print(barra1.fep)
#rint(barra2.fep)


#barra1.comprimento_barra()
#barra1.calcula_klocal()
#barra1.calcula_r()
#barra1.calcula_k()
#barra1.calcula_q()


#print(barra1.kl)
#print(barra1.r)
#print(barra1.k)
#print(barra1.q)


#carregamento1.mostre()
#carregamento1.mostre()
#carregamento1.mostre()

#    E = 10000.0
#    area = 2.0 * 5.0
#    inercia = 1000.0











    





    




