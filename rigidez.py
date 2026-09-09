#import matplotlib.pyplot as plt
import numpy as np
import math
#from matplotlib.collections import LineCollection

class Concreto:
    """Representa as propriedades e parâmetros referente ao concreto."""

    def __init__(self, fck: float, γc: float = 1.4, αE: float = 1.0) -> None:
        """
        Instancia um objeto da classe Concreto


        Args:

            fck: resistência característica do concreto
            γc: coeficiente de segurança do concreto
            αE: coeficiente de ajuste do módulo de elasticidade

            Valores de αE:

            αE = 1.2 para basalto e diabásio
            αE = 1.0 para granito e gnaisse
            αE = 0.9 para calcário
            αE = 0.7 para arenito

        Raises:
            ValueError: Se o fck estiver fora do intervalo de 20 a 90 MPa.
            ValueError: Caso αE não corresponda aos valores indicados.
        """

        if fck > 90 or fck < 20:
            raise ValueError("fck deve ser entre 20 e 90")
        
        if αE not in [0.7, 0.9, 1.0, 1.2]:
            raise ValueError("αE deve ser 0.7, 0.9, 1.0 ou 1.2")
        
        self.fck = fck 
        self.γc = γc
        self.fcd = fck/γc
        self.Eci = αE*5600*math.sqrt(fck)
        self.αc = 0.85
        self.λ = 0.8
        self.εcu = 0.0035
        self.εc2 = 0.002

    @property
    def ηc(self) -> float:
        ηc = 1.0
        if self.fck > 40 and self.fck <= 90:
            ηc = (40/self.fck)**(1/3)
        return ηc

    @property
    def qualquer(self):
        if self.fck > 50 and self.fck <= 90:
            self.αc = 0.85*(1-(self.fck-50)/200)
            self.εcu = 0.0026 + 0.035*((90-self.fck)/100)**4
            self.εc2 = 0.002 + 0.000085*((self.fck-50)**0.53)
            self.Eci = 21.5*10**3*self.αE*((self.fck/10)+1.25)**(1/3)
            self.λ = 0.8 - ((self.fck-50)/400)


class Aco:
    def __init__(self, fyk:float, γs:float = 1.15, Es:float = 210000.0):
    
        """
        Cria o objeto Aço
    
        Argumentos:
    
            fyk: resistência característica do aço
            γs: coeficiente de segurança do aço
            Es: módulo de elasticidade do aço
        """
    
        self.fyk = fyk 
        self.γs = γs
        self.fyd = fyk/γs
        self.Es = Es
        self.εyd = self.fyd/Es


class Secao:
    def __init__(self,nome:str, b:float, h:float, E:float = None, d:float = None, d_linha:float = None):
    
        """
        Cria o objeto Seção
    
        Argumentos:
    
            nome: nome da seção
            b: largura da seção
            h: altura da seção
            E: módulo de elasticidade da seção
            d: altura útil da armadura tracionada
            d_linha: altura útil da armadura comprimida     
        """
    
        self.nome = nome
        self.b = b
        self.h = h
        self.d = d
        self.d_linha = d_linha
        self.E = E
        self.A = b*h
        self.I = b*h**3/12

class No:
    def __init__(self, num:int, x:float, y:float):

        """
        Cria o objeto no

        Argumentos:

            num: numero do no
            x: coordenada x do no
            y: coordenada y do no     
        """

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
    def __init__(self, num:int, no1:No, no2:No, secao:Secao, concreto:Concreto = None, aco:Aco = None):

        """
        Cria o objeto barra

        Argumentos:

            num: numero da barra
            no1: no inicial
            no2: no final
            concreto: objeto com propriedades do concreto
            aco: objeto com propriedades do aço
            secao: objeto com propriedades da seção transversal      
        """

        self.num = int(num)
        self.noi = no1
        self.noj = no2
        self.concreto = concreto
        self.aco = aco
        self.secao = secao
        self.L = None      
        self.kl = np.zeros((6,6))
        self.k = np.zeros((6,6))
        self.r = np.zeros((6,6))
        self.fepl = np.zeros((6))
        self.fep = np.zeros((6))
        self.f = np.zeros((6))
        self.fl = np.zeros((6))
        self.u = np.zeros((6))
        self.q = np.zeros((6),dtype=int)
        self.x = np.zeros((2))
        self.ASL = np.zeros((2))
        self.comprimento_barra()
        self.calcula_r()
        
    def comprimento_barra(self) -> float:
        dx = self.noj.x - self.noi.x
        dy = self.noj.y - self.noi.y
        self.L = math.sqrt( dx * dx + dy * dy )

    def calcula_klocal(self) -> float:
        """
        Calcula matriz de rigidez no sistema local para a barra      
        """
        
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

    def calcula_r(self) -> float:
        """
        Calcula matriz de rotação para a barra      
        """ 

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

    def calcula_k(self) -> float:
        """
        Monta matriz de rigidez no sistema global para a barra      
        """

        self.k = np.linalg.inv(self.r)@ self.kl @ self.r

    def monta_q(self) -> int:
        """
        Monta vetor de correspondencia de graus de liberdade da barra      
        """

        z = -1
        M = np.array([self.noi.num, self.noj.num])
        for j in range(0,2):
            for jk in range(0,3):
                z = z + 1
                self.q[z] = 3*(M[j]-1)+jk

    def calcula_vetor_ASL(self) -> int:
        """
        Calcula vetor x e  vetor ASL da barra      
        """

        γq = 1.4
        M = np.array([self.fl[2], self.fl[5]])
        N = np.array([self.fl[0], self.fl[3]])
        Md = M*γq
        Nd = N*γq

        xlim = (self.concreto.εcu*self.secao.d)/(self.concreto.εyd+self.concreto.εcu)

        if xlim > 0.45*self.secao.d:
            xlim = 0.45*self.secao.d
        if self.concreto.fck > 50 and self.concreto.fck <= 90:
            if xlim > 0.35*self.secao.d:
                xlim = 0.35*self.secao.d

        for i in N:
            if N[i] == 0:
                Mdlim = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*xlim*self.concreto.λ(self.secao.d-self.concreto.λ*xlim/2.0)
                if Md[i] > Mdlim:
                    self.x[i] = (self.secao.d - math.sqrt(self.secao.d**2-2.0*Md[i]/(self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b)))/self.concreto.λ
                    self.ASL[i] = (self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.x[i]*self.concreto.λ)/self.aco.fyd



    








            

        
      
class Estrutura:
    def __init__(self, nos:list, barras:list):

        """
        Cria o objeto estrutura

        Argumentos:

            nos: lista de objetos nós
            barras: lista de objetos barras      
        """

        self.nos = nos
        self.barras = barras
        self.nnos = len(nos)
        self.fnos = np.zeros((self.nnos*3))
        self.f = np.zeros((self.nnos*3))
        self.u = np.zeros((self.nnos*3))
        self.k = np.zeros((self.nnos*3,self.nnos*3))
        self.k01 = np.zeros((self.nnos*3,self.nnos*3))
        self.R = np.zeros((self.nnos*3))
              
    def monta_k(self) -> float:
        """
        Monta a matriz de rigidez da estrutura       
        """

        for barra in self.barras:
            barra.comprimento_barra()
            barra.calcula_klocal()
            barra.calcula_r()
            barra.calcula_k()
            barra.monta_q()
            for j in range(0,6):
                for jk in range(0,6):
                    self.k[barra.q[j], barra.q[jk]] += barra.k[j,jk]

    def monta_k01(self) -> float:
        """
        Aplica as condições de contorno na matriz da estrutura para calcular os deslocamentos       
        """

        self.k01 = self.k
        for no in self.nos:
            if no.Tx:
                gdl = 3*no.num-3
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0
            if no.Ty:
                gdl = 3*no.num-2
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0               
            if no.Rz:
                gdl = 3*no.num-1
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0           
        
    def monta_fnos(self) -> float:
        """
        Monta vetor de forças necessário para calcular deslocamentos na estrutura      
        """

        for no in self.nos:
            self.fnos[3*no.num-3] = no.Fx
            self.fnos[3*no.num-2] = no.Fy
            self.fnos[3*no.num-1] = no.Mz 
            
        for barra in self.barras:
            for j in range(0,6):
                self.fnos[barra.q[j]] += -barra.fep[j]
                
    def aplica_cc_fnos(self) -> float:
        """
        Aplica as condições de contorno no vetor fnos para calcular deslocamentos na estrutura       
        """
        for no in self.nos:
            if no.Tx:
                self.fnos[3*no.num-3] = 0
            if no.Ty:
                self.fnos[3*no.num-2] = 0
            if no.Rz:
                self.fnos[3*no.num-1] = 0

    def calcula_deslocamentos(self) -> float:
        """
        Calcula o vetor de deslocamentos da estrutura       
        """

        self.u = np.linalg.inv(self.k01)@ self.fnos

    def calcula_solicitacoes_internas_nodais(self) -> float:
        """
        Calcula as forças no sistema local para as barras = solicitações
        """

        for barra in self.barras:
            barra.monta_q()
            barra.calcula_r()

            """
            Monta o vetor u para as barras
            """
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]

            """
            Calcula o vetor f para as barras
            """      
            barra.f = barra.k @ barra.u + barra.fep


            """
            Transforma o vetor f para o sistema local para as barras = solicitacoes
            """       
            barra.fl = barra.r @ barra.f    

    def calcula_reacoes(self) -> float:
        """
        Calcula um vetor com as reações da estrutura       
        """

        for barra in self.barras:
            barra.monta_q()

            """
            Monta o vetor u para as barras
            """
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]
            
            """
            Calcula o vetor f para as barras       
            """
            barra.f = barra.k @ barra.u + barra.fep

            if barra.noi.Tx:
                self.R[3*barra.noi.num-3] += barra.f[0]
            if barra.noi.Ty:
                self.R[3*barra.noi.num-2] += barra.f[1]
            if barra.noi.Rz:
                self.R[3*barra.noi.num-1] += barra.f[2]

            if barra.noj.Tx:
                self.R[3*barra.noj.num-3] += barra.f[3]
            if barra.noj.Ty:
                self.R[3*barra.noj.num-2] += barra.f[4]
            if barra.noj.Rz:
                self.R[3*barra.noj.num-1] += barra.f[5]

        for no in self.nos:

            if no.Tx:
                self.R[3*no.num-3] += -no.Fx
            if no.Ty:
                self.R[3*no.num-2] += -no.Fy
            if no.Rz:
                self.R[3*no.num-1] += -no.Mz

class Carregamento_distribuido:
     
    def __init__(self, a:float, lw:float, w1:float, w2:float, barra:Barra):

        """
        Cria o objeto carregamento distribuido

        Argumentos:

            a: distância do nó inicial da barra até o inicio do carregamento
            lw: comprimento do carregamento
            w1: modulo inicial do carregamento
            w2: modulo final do carregamento
            barra: objeto barra
        """

        self.a = a
        self.lw = lw
        self.w1 = w1
        self.w2 = w2
        self.barra = barra
        
    def calcula_fepl(self) -> float:
        """
        Calcula o vetor de forças de engastamento perfeito no sistema local da barra
        """

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

    def calcula_fep(self) -> float:
        """
        Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra
        """

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl

class Carregamento_pontual:
   
    def __init__(self, a:float, Px:float, Py:float, barra:Barra):
        
        """
        Cria o objeto carregamento pontual

        Argumentos:

            a: distância do nó inicial da barra até o ponto de aplicação da força
            Px: componente x da força
            Py: componente y da força
            barra: objeto barra
        """

        self.a = a
        self.Px = Px
        self.Py = Py
        self.barra = barra

    def calcula_fepl(self) -> float:
        """
        Calcula o vetor de forças de engastamento perfeito no sistema local da barra
        """

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

    def calcula_fep(self) -> float:
        """
        Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra
        """

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl

class Modelo:

    def __init__(self, bw:float, h:float, l:float, lf:float, a:float, b:float, c:float, ap1:float, ap2:float, P:float, w:float):

        """
        Inicializa o modelo estrutural

        Argumentos:

            bw: largura viga
            h: altura da viga
            l: comprimento da viga
            lf: comprimento até centro da abertura
            a: altura da abertura
            b: largura da abertura
            c: espessura do banzo inferior na regiao da abertura
            ap1: largura pilar esquerdo
            ap2: largura pilar direito
            P: carregamento pontual
            w: modulo carregamento distribuido
        """

        self.bw = bw
        self.h = h
        self.l = l
        self.lf = lf
        self.a = a
        self.b = b
        self.c = c
        self.ap1 = ap1
        self.ap2 = ap2
        self.P = P
        self.w = w
        self.nos = None
        self.barras = None

        """
        Definição dos nós
        """
        no1 = No(1, 0.0, 0.0)
        no2 = No(2, self.ap1/2.0, 0.0)
        no3 = No(3, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, 0.0)
        no4 = No(4, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, -(self.h-self.c)/2.0)
        no5 = No(5, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, (self.a+self.c)/2.0)
        no6 = No(6, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, 0.0)
        no7 = No(7, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, -(self.h-self.c)/2.0)
        no8 = No(8, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, (self.a+self.c)/2.0)
        no9 = No(9, self.ap1/2.0 + self.l, 0.0)
        no10 = No(10, self.ap1/2.0 + self.l + self.ap2/2.0, 0.0)
        self.nos = [no1, no2, no3, no4, no5, no6, no7, no8, no9, no10]

        """
        Aplicação das restrições nodais
        """
        no1.Tx = True
        no1.Ty = True
        no1.Rz = True
        no10.Tx = True
        no10.Ty = True
        no10.Rz = True

        """
        Definição das barras e propriedades
        """
        E = 280000.0
        I1 = self.bw*self.h**3/12.0
        A1 = self.bw*self.h

        I2 = self.bw*self.c**3/12.0
        A2 = self.bw*self.c

        I3 = self.bw*(self.h-self.a-self.c)**3/12.0
        A3 = self.bw*(self.h-self.a-self.c)

        I4 = I1*1000
        A4 = A1*1000

        """
        Ordem dos nós na definição das barras de ser sempre da esquerda para direita ou de baixo para cima
        """
        barra1 = Barra(1,no1,no2,E,A1,I1)
        barra2 = Barra(2,no2,no3,E,A1,I1)
        barra3 = Barra(3,no4,no7,E,A2,I2)
        barra4 = Barra(4,no5,no8,E,A3,I3)
        barra5 = Barra(5,no6,no9,E,A1,I1)
        barra6 = Barra(6,no9,no10,E,A1,I1)
        barra7 = Barra(7,no4,no3,E,A4,I4)
        barra8 = Barra(8,no3,no5,E,A4,I4)
        barra9 = Barra(9,no7,no6,E,A4,I4)
        barra10 = Barra(10,no6,no8,E,A4,I4) 
        self.barras = [barra1, barra2, barra3, barra4, barra5, barra6, barra7, barra8, barra9, barra10]

        """
        Definição dos carregamentos nas barras
        """
        carregamento1 = Carregamento_distribuido(0.0,self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra2)
        carregamento2 = Carregamento_distribuido(0.0,self.b+self.bw,self.w,self.w,barra4)
        carregamento3 = Carregamento_distribuido(0.0,self.l-self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra5)
        
        """
        Calcula as forças de engastamento perfeito
        """
        carregamento1.calcula_fepl()
        carregamento1.calcula_fep()
        carregamento2.calcula_fepl()
        carregamento2.calcula_fep()
        carregamento3.calcula_fepl()
        carregamento3.calcula_fep()
        
        """
        Define a estrutura, calcula os resultados de deslocamentos, solicitações internas e reações
        """
        portico = Estrutura(self.nos,self.barras)
        portico.monta_k()
        print("k")
        print_matriz(portico.k)

        portico.monta_k01()
        print("k01")
        print_matriz(portico.k01)

        portico.monta_fnos()
        print("fnos")
        print_vetor(portico.fnos)

        portico.aplica_cc_fnos()
        print("f01")
        print_vetor(portico.fnos)

        portico.calcula_deslocamentos()
        print("deslocamentos")
        print_vetor(portico.u)

        portico.calcula_solicitacoes_internas_nodais()
        print("solicitacoes barra1")
        print_vetor(portico.barras[0].fl)
        print("solicitacoes barra2")
        print_vetor(portico.barras[1].fl)
        print("solicitacoes barra3")
        print_vetor(portico.barras[2].fl)
        print("solicitacoes barra4")
        print_vetor(portico.barras[3].fl)
        print("solicitacoes barra5")
        print_vetor(portico.barras[4].fl)
        print("solicitacoes barra6")
        print_vetor(portico.barras[5].fl)
        print("solicitacoes barra7")
        print_vetor(portico.barras[6].fl)
        print("solicitacoes barra8")
        print_vetor(portico.barras[7].fl)
        print("solicitacoes barra9")
        print_vetor(portico.barras[8].fl)
        print("solicitacoes barra10")
        print_vetor(portico.barras[9].fl)

        portico.calcula_reacoes()
        print("reacoes")
        print_vetor(portico.R)

def print_matriz(m):
    l, c = m.shape

    tamanho_colunas = []
    for i in range(c):
        maior_l = 0
        for j in range(l):
            s = f"{m[j,i]:.1f}"
            maior_l = max(len(s), maior_l)
        tamanho_colunas.append(maior_l)

    cabecalho = '\x1b[90m  '
    for i in range(c):
        cabecalho += ' ' + str(i).center(tamanho_colunas[i])
    cabecalho += '\x1b[0m'
    print(cabecalho)

    for i in range(l):
        linha = '\x1b[90m' + str(i).rjust(2) + '\x1b[0m'
        for j in range(c):
            s = f"{m[j,i]:.1f}".rjust(tamanho_colunas[j])
            linha += ' ' + s
        print(linha)


def print_vetor(v):
    def fmt(n):
        if abs(n) < 0.06 or abs(n) > 1e6:
            return f'{n:#.2g}'
        return f'{n:.1f}'

    l = v.shape[0]

    maior_l = 0
    for i in range(l):
        s = fmt(v[i])
        maior_l = max(len(s), maior_l)

    for i in range(l):
        linha = '\x1b[90m' + str(i).rjust(2) + '\x1b[0m'
        s = fmt(v[i]).rjust(maior_l)
        linha += ' ' + s
        print(linha)












    





    




