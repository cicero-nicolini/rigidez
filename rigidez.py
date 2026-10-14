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

            fck: resistência característica do concreto [MPa]
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
        self.αE = αE
    
    @property
    def fct_m(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fct_m

    @property
    def fctk_sup(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_sup = 1.3*fct_m
        return fctk_sup
        
    @property
    def fctk_inf(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_inf = 0.7*fct_m
        return fctk_inf

    @property
    def fctd(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_inf = 0.7*fct_m
        fctd = fctk_inf/self.γc
        return fctd

    @property
    def αv2(self) -> float:
        αv2 = (1 - (self.fck/250))
        return αv2 

    @property
    def fcd(self) -> float:
        fcd = (self.fck/self.γc)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fcd
    
    @property
    def ηc(self) -> float:
        ηc = 1.0
        if 40 < self.fck <= 90:
            ηc = (40/self.fck)**(1/3)
        return ηc
    
    @property
    def Eci(self) -> float:
        Eci = (5600*self.αE*math.sqrt(self.fck))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            Eci = (21.5*10**3*self.αE*((self.fck/10)+1.25)**(1/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return Eci

    @property
    def αc(self) -> float:
        αc = 0.85
        if 50 < self.fck <= 90:
            αc = 0.85*(1-(self.fck-50)/200)
        return αc

    @property
    def λ(self) -> float:
        λ = 0.8
        if 50 < self.fck <= 90:
            λ = 0.8 - ((self.fck-50)/400)
        return λ

    @property
    def εc2(self) -> float:
        εc2 = 0.002
        if 50 < self.fck <= 90:
            εc2 = 0.002 + 0.000085*((self.fck-50)**0.53)
        return εc2

    @property
    def εcu(self) -> float:
        εcu = 0.0035
        if 50 < self.fck <= 90:
            εcu = 0.0026 + 0.035*((90-self.fck)/100)**4
        return εcu
  
class Aco:
    """Representa as propriedades e parâmetros referente ao aço."""

    def __init__(self, fyk:float, fywk:float, γs:float = 1.15, Es:float = 21000.0) -> None:
        """
        Instancia um objeto da classe Aco


        Args:

            fyk: resistência característica do aço [MPa]
            γs: coeficiente de segurança do aço
            Es: módulo de elasticidade do aço
        """

        self.fyk = fyk
        self.fywk = fywk 
        self.γs = γs
        self.Es = Es #Unidade compatibilizada para testes de MPa para kN/cm²
    
    @property
    def fyd(self) -> float:
        fyd = (self.fyk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fyd

    @property
    def fywd(self) -> float:
        fywd = (self.fywk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fywd
    
    @property
    def εyd(self) -> float:
        fyd = (self.fyk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        εyd = fyd/self.Es
        return εyd

class Secao:
    """Representa as propriedades e parâmetros referente a seção transversal."""

    def __init__(self,nome:str, b:float, h:float, E:float = None, d:float = None, d_linha:float = None) -> None:
        """
        Instancia um objeto da classe Secao


        Args:

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
        self.yt = h/2 #distância do centro de gravidade até a fibra mais tracionada
        self.W0 = (b*h**3/12)/(h/2) #módulo de resistência à flexão da seção bruta de concreto 
                                    #em relação à fibra mais tracionada

class No:
    """Representa as propriedades e parâmetros referente ao objeto nó."""

    def __init__(self, num:int, x:float, y:float) -> None:
        """
        Instancia um objeto da classe Nó


        Args:

            num: número do nó
            x: coordenada x do nó
            y: coordenada y do nó
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
    """Representa as propriedades e parâmetros referente ao objeto barra."""

    def __init__(self, num:int, no1:No, no2:No, secao:Secao, concreto:Concreto = None, aco:Aco = None, modelo:int = 1, α:float = 90.0, θ:float = 30.0) -> None:
        """
        Instancia um objeto da classe Barra


        Args:

            num: numero da barra
            no1: no inicial
            no2: no final
            secao: objeto com propriedades da seção transversal 
            concreto: objeto com propriedades do concreto
            aco: objeto com propriedades do aço
            modelo: modelo de cálculo da armadura transversal (modelo 1 ou modelo 2)
            α: inclinação dos estribos (α deve estar entre 45° e 90°)
            θ: inclinação da biela comprimida de concreto (necessário para modelo 2, θ deve estar entre 30° e 45°)

        Raises:
            ValueError: Caso o modelo não corresponda aos valores indicados.
            ValueError: Se α estiver fora do intervalo indicado.
            ValueError: Se θ estiver fora do intervalo indicado.
        """

        if modelo not in [1, 2]:
            raise ValueError("modelo deve ser 1 ou 2")
        if not 45 <= α <= 90:
            raise ValueError("Erro: α deve estar entre 45° e 90°.")
        if not 30 <= θ <= 45:
            raise ValueError("Erro: θ deve estar entre 30° e 45°.")

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
        self.ASL2 = np.zeros((2))
        self.ASW_S = np.zeros((2))
        self.comprimento_barra()
        self.calcula_r()
        self.modelo = modelo
        self.α = α
        self.θ = θ
        
    def comprimento_barra(self) -> None:
        """Calcula o comprimento da barra."""

        dx = self.noj.x - self.noi.x
        dy = self.noj.y - self.noi.y
        self.L = math.sqrt( dx * dx + dy * dy )

    def calcula_klocal(self) -> None:
        """Calcula a rigidez da barra no sistema local da barra."""

        # Calculando os termos da matriz de rigidez  
        a1 = self.secao.E*self.secao.A/self.L
        a2 = (12.0*self.secao.E*self.secao.I)/self.L**3
        a3 = (6.0*self.secao.E*self.secao.I)/self.L**2
        a4 = (4.0*self.secao.E*self.secao.I)/self.L
        a5 = (2.0*self.secao.E*self.secao.I)/self.L

        # Preenchendo a matriz de rigidez local
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

    def calcula_r(self) -> None:
        """Calcula matriz de rotação para a barra."""        

        # Calculando termos da matriz de rotação
        dx = self.noj.x - self.noi.x
        dy = self.noj.y - self.noi.y
        s = dy/self.L
        c = dx/self.L

        # Preenchendo a matriz de rotação
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

    def calcula_k(self) -> None:
        """Monta matriz de rigidez no sistema global para a barra."""          

        self.k = np.linalg.inv(self.r)@ self.kl @ self.r

    def monta_q(self) -> None:
        """Monta a função vetor q de correspondencia da barra."""       

        # Matriz que guarda a numeração dos nós
        M = np.array([self.noi.num, self.noj.num])

        # Monta o vetor de correspondencia q de acordo com a referncia
        # [Soriano, 2005, p. 76] considerando g = 3
        z = -1                  # z = -1 para ajustar a contagem
        for j in range(0,2):
            for jk in range(0,3):
                z = z + 1
                self.q[z] = 3*(M[j]-1)+jk

    def calcula_vetor_ASW_S(self, modelo:int, α, θ, M, N, V) -> None:
        """
        Calcula vetor ASW_S que resulta na área de armadura transversal em cm²/m.
        
        Raises
            ValueError: Se a força cortante solicitante de cálculo (Vsd) for superior
            à resistência das bielas comprimidas de concreto (Vrd2).
        """

        α = math.radians(α)
        θ = math.radians(θ)

        # Definindo o coeficiente majorador de esforços
        γf = 1.4

        # Guardando as solicitações (M, N) dos nós extremos da barra
        #M = np.array([self.fl[2], self.fl[5]])
        #N = np.array([self.fl[0], self.fl[3]])
        #V = np.array([self.fl[1], self.fl[4]])

        Vsd = abs(V*γf)

        for i in range(len(V)):

            if modelo == 1:
                Vrd2 = 0.27*self.concreto.αv2*self.concreto.fcd*self.secao.b*self.secao.d
    
                if Vsd[i] <= Vrd2:
                    Vc0 = 0.6*self.concreto.fctd*self.secao.b*self.secao.d

                    if N[i] >= 0 and self.x[i] < 0 or N[i] >= 0 and self.x[i] > self.secao.h:
                        Vc = 0
                    else:
                        Vc = Vc0
                    
                    Vsw = Vsd[i] - Vc
                    self.ASW_S[i] = (Vsw/(0.9*self.secao.d*self.aco.fywd*(math.sin(α) + math.cos(α))))*100 #Compatibilização de unidades para testes de cm²/cm para cm²/m 
                    ASW_S_min = ((0.2*self.concreto.fct_m*self.secao.b*math.sin(α))/self.aco.fywk)*1000 #Compatibilização de unidades para testes para cm²/m 
                    self.ASW_S[i] = max(self.ASW_S[i], ASW_S_min)
                else:
                    raise ValueError("Vsd é maior que Vrd2, revise as dimensões da seção ou o concreto.")

            if modelo == 2:
                Vrd2 = 0.54*self.concreto.αv2*self.concreto.fcd*self.secao.b*self.secao.d*math.sin(θ)**2*((1/math.tan(θ)) + (1/math.tan(α)))

                if Vsd[i] <= Vrd2:
                    Vc0 = (0.42*self.concreto.fct_m*self.secao.b*self.secao.d)/self.concreto.γc
                    Vc1 = ((Vrd2 - Vsd[i])/(Vrd2 - Vc0))*Vc0

                    if N > 0 and self.x[i] > self.secao.h:
                        Vc = 0
                    else:
                        Vc = Vc1

                    self.ASW_S[i] = ((Vsd[i] - Vc)/(0.9*self.secao.d*self.aco.fywd*math.sin(α)*((1/math.tan(θ)) + (1/math.tan(α)))))*100 #Compatibilização de unidades para testes de cm²/cm para cm²/m 
                    ASW_S_min = ((0.2*self.concreto.fct_m*self.secao.b*math.sin(α))/self.aco.fywk)*1000 #Compatibilização de unidades para testes para cm²/m 
                    self.ASW_S[i] = max(self.ASW_S[i], ASW_S_min)
                else:
                    raise ValueError("Vsd é maior que Vrd2, revise as dimensões da seção ou o concreto.")
                
    def calcula_vetor_ASL(self) -> None:
        """Calcula vetor x, vetor ASL e ASL2 para a barra."""

        # Definindo o coeficiente majorador de esforços
        γf = 1.4

        # Guardando as solicitações (M, N) dos nós extremos da barra
        M = np.array([self.fl[2], self.fl[5]])
        N = np.array([self.fl[0], self.fl[3]])

        # Majorando os esforços
        Md = abs(M*γf)
        Nd = abs(N*γf)

        # Calculando o xlim
        xlim = (self.concreto.εcu*self.secao.d)/(self.aco.εyd+self.concreto.εcu)

        # Verificando o xlim de acordo com o limite de ductilidade conforme  
        # a referencia [NBR6118:2026, p. x]
        if xlim > 0.45*self.secao.d:
            xlim = 0.45*self.secao.d

        if self.concreto.fck > 50 and self.concreto.fck <= 90:
            if xlim > 0.35*self.secao.d:
                xlim = 0.35*self.secao.d

        # Calculando o Mdlim
        Mdlim = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*xlim*self.concreto.λ*(self.secao.d - self.concreto.λ*xlim/2.0)

        Md_min = 0.8*self.secao.W0*self.concreto.fctk_sup
        xmin = (self.secao.d - math.sqrt(self.secao.d**2-2.0*Md_min/(self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b)))/self.concreto.λ
        Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*xmin
        ASL_min = Rcc/self.aco.fyd
        
        # Para cada i e n, onde i é o índice do vetor N e n é o valor
        # correspondente ao esforço normal armazenado na posição i do vetor N.
        for i, n in enumerate(N):
            
            # Verificação do tipo de flexão: simples ou composta
            if n == 0:
                # Flexão simples
                
                if Md[i] < Mdlim:
                    # Domínio 2 ou 3
                    
                    self.x[i] = (self.secao.d - math.sqrt(self.secao.d**2-2.0*Md[i]/(self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b)))/self.concreto.λ
                    Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]
                    self.ASL[i] = Rcc/self.aco.fyd
                    self.ASL[i] = max(self.ASL[i], ASL_min)
                else:
                    # armadura dupla
                   
                    # Limite do domínio 3
                    self.x[i] = xlim
                    Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]
                    z = self.secao.d - (self.concreto.λ*self.x[i]/2.0)
                    ε2 = self.concreto.εcu*(1-(self.secao.d_linha/self.x[i]))

                    # Verificando o regime da armadura de compressão
                    if ε2 > self.aco.εyd:
                        σ2 = self.aco.fyd
                    else:
                        σ2 = self.aco.Es*ε2
                    # Calculo da área de aço de compressão e tração
                    self.ASL2[i] = (Md[i] - Rcc*z)/(σ2*(self.secao.d-self.secao.d_linha))
                    self.ASL[i] = (Rcc + self.ASL2[i]*σ2)/self.aco.fyd
                    self.ASL[i] = max(self.ASL[i], ASL_min)

            if n > 0:
                # Flexo-tração
                e0 = M/N  # calculo da excentricidade da carga
                
                # Verificando regime da flexo-tração
                if e0[i] >= (self.secao.d-self.secao.d_linha)/2.0:
                    # Flexo-tração com grande excentricidade (domínio 2 ou 3)
                    
                    # Calculo da excentricidade da armadura 1
                    e1 = e0[i] - (self.secao.d-self.secao.d_linha)/2.0

                    # Verificando o tipo de armadura
                    if Nd[i]*e1 <= Mdlim:
                        # armadura simples
                        
                        # Resolvendo o sistema de equilíbrio para encontrar o x
                        a = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ**2/2.0
                        b = -self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.d*self.concreto.λ
                        c = Nd[i]*e1
                        delta = b**2-4*a*c
                        x1 = (-b+math.sqrt(delta))/(2*a)
                        x2 = (-b-math.sqrt(delta))/(2*a)
                        if x1 >= 0 and x1 <= xlim:
                            self.x[i] = x1
                        if x2 >= 0 and x2 <= xlim:
                            self.x[i] = x2

                        Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]

                        # Calculo da área de aço armadura 1
                        self.ASL[i] = (Nd[i] + Rcc)/self.aco.fyd
                        self.ASL[i] = max(self.ASL[i], ASL_min)
                    else:
                        # armadura dupla

                        self.x[i] = xlim
                        Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]
                        
                        # Calculo da tensão na armadura 2
                        ε2 = self.concreto.εcu*(self.x[i]-self.secao.d_linha/self.x[i])
                    
                        if ε2 > self.aco.εyd:
                            σ2 = self.aco.fyd
                        else:
                            σ2 = self.aco.Es*ε2

                        # Calculo da armadura 2
                        self.ASL2[i] = (Nd[i]*e1 - Mdlim)/(σ2*(self.secao.d-self.secao.d_linha))
                        self.ASL[i] = (Nd[i] + self.ASL2[i]*σ2 + Rcc)/(self.aco.fyd)
                        self.ASL[i] = max(self.ASL[i], ASL_min)
                else:
                    # Flexo-tração com pequena excentricidade (domínio 1)
                    
                    e1 = (self.secao.d-self.secao.d_linha)/2.0 - e0[i]
                    e2 = (self.secao.d-self.secao.d_linha)/2.0 + e0[i]
                    σ1 = self.aco.fyd
                    σ2 = σ1

                    self.ASL2[i] = (Nd[i]*e1)/(σ2*(self.secao.d-self.secao.d_linha))
                    self.ASL[i] = (Nd[i]*e2)/(σ1*(self.secao.d-self.secao.d_linha))
                    self.ASL[i] = max(self.ASL[i], ASL_min)
                    
            if n < 0:
                # Flexo-compressão
               
                e0 = abs(M/N)  # calculo da excentricidade da carga
                
                e1 = (self.secao.d-self.secao.d_linha)/2.0 + e0[i]
                e2 = (self.secao.d-self.secao.d_linha)/2.0 - e0[i]
                #verificação da necessidade teórica de armadura
                e2_0 = Nd[i]/(1.7*self.concreto.fcd*self.secao.b) - self.secao.d_linha
                e2_gp = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*xlim*((self.concreto.λ*xlim/2.0)-self.secao.d_linha)/Nd[i]

                if e2_0 <= e2:
                    # Não há necessidade de armadura
                    # ASL_min
                    
                    self.ASL[i] = ASL_min
                else:
                    # Necessidade de armadura
                       
                    if e2 < 0 or e2 < e2_gp:
                        # Flexo-compressão com grande excentricidade
                        
                        self.x[i] = xlim
                       
                        if Nd[i]*e1 <= Mdlim:
                            # armadura simples
                            
                            a = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ**2/2.0
                            b = -self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.d*self.concreto.λ
                            c = Nd[i]*e1
                            
                            delta = b**2-4*a*c
                            x1 = (-b+math.sqrt(delta))/(2*a)
                            x2 = (-b-math.sqrt(delta))/(2*a)
                            if x1 >= 0 and x1 <= xlim:
                                self.x[i] = x1
                            if x2 >= 0 and x2 <= xlim:
                                self.x[i] = x2
                            
                            Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]
                            self.ASL[i] = (Rcc-Nd[i])/self.aco.fyd
                            self.ASL[i] = max(self.ASL[i], ASL_min)
                        else:
                            # armadura dupla
                            
                            self.x[i] = xlim
                            ε2 = self.concreto.εcu*(xlim-self.secao.d_linha/xlim)
                            if ε2 < self.aco.εyd:
                                σ2 = self.aco.Es*ε2
                            else:
                                σ2 = self.aco.fyd
                            self.ASL2[i] = (Nd[i]*e1 - Mdlim)/(σ2*(self.secao.d-self.secao.d_linha))
                            self.ASL[i] = (self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*xlim - Nd[i] + self.ASL2[i]*σ2)/self.aco.fyd
                            self.ASL[i] = max(self.ASL[i], ASL_min)
                    else:
                        e2_pc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.h*(0.5*self.secao.h-self.secao.d_linha)/Nd[i]
                        
                        if e2 < e2_pc:
                            # Flexo-compressão com pequena excentricidade
                            
                            a = -self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ**2/2.0
                            b = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.d_linha*self.concreto.λ
                            c = Nd[i]*e2
                            delta = b**2-4*a*c
                            x1 = (-b+math.sqrt(delta))/(2*a)
                            x2 = (-b-math.sqrt(delta))/(2*a)
                            if x1 >  xlim:
                                self.x[i] = x1
                            if x2 >  xlim:
                                self.x[i] = x2

                            Rcc = self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.concreto.λ*self.x[i]

                            if xlim < self.x[i] <= self.secao.h:
                                ε2 = self.concreto.εcu*(self.x[i]-self.secao.d_linha/self.x[i])
                            else:
                                ε2 = (2/1000)*(self.x[i]-self.secao.d_linha)/(self.x[i]-3*self.secao.h/7)
                            if ε2 < self.aco.εyd:
                                σ2 = self.aco.Es*ε2
                            else:
                                σ2 = self.aco.fyd
                            self.ASL2[i] = (Nd[i] - Rcc)/σ2
                            self.ASL[i] = ASL_min
                        else:
                            # compressão composta
                            
                            ε1 = 2/1000
                            ε2 = ε1
                            if ε2 <= self.aco.εyd:
                                σ2 = self.aco.Es*ε2
                            else:
                                σ2 = self.aco.fyd
                            σ1 = σ2
                            self.ASL2[i] = (Nd[i]*e1 - self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.h*(self.secao.d - self.secao.h*0.5))/(σ2*(self.secao.d-self.secao.d_linha))
                            self.ASL[i] = (Nd[i]*e2 - self.concreto.αc*self.concreto.ηc*self.concreto.fcd*self.secao.b*self.secao.h*(self.secao.h*0.5-self.secao.d_linha))/(σ1*(self.secao.d-self.secao.d_linha))
                            self.ASL[i] = max(self.ASL[i], ASL_min)
                                                         
class Estrutura:
    """Representa as propriedades e parâmetros referente ao objeto estrutura."""

    def __init__(self, nos:list, barras:list) -> None:
        """
        Instancia um objeto da classe Estrutura


        Args:

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
              
    def monta_k(self) -> None:
        """Monta a matriz de rigidez da estrutura."""        
        
        for barra in self.barras:
            barra.comprimento_barra()
            barra.calcula_klocal()
            barra.calcula_r()
            barra.calcula_k()
            barra.monta_q()
            for j in range(0,6):
                for jk in range(0,6):
                    self.k[barra.q[j], barra.q[jk]] += barra.k[j,jk]

    def monta_k01(self) -> None:
        """Aplica as condições de contorno na matriz da estrutura para calcular os deslocamentos."""     
    
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
        
    def monta_fnos(self) -> None:
        """Monta vetor de forças necessário para calcular deslocamentos na estrutura."""    

        for no in self.nos:
            self.fnos[3*no.num-3] = no.Fx
            self.fnos[3*no.num-2] = no.Fy
            self.fnos[3*no.num-1] = no.Mz 
            
        for barra in self.barras:
            for j in range(0,6):
                self.fnos[barra.q[j]] += -barra.fep[j]
                
    def aplica_cc_fnos(self) -> None:
        """Aplica as condições de contorno no vetor fnos para calcular deslocamentos na estrutura."""      
    
        for no in self.nos:
            if no.Tx:
                self.fnos[3*no.num-3] = 0
            if no.Ty:
                self.fnos[3*no.num-2] = 0
            if no.Rz:
                self.fnos[3*no.num-1] = 0

    def calcula_deslocamentos(self) -> None:
        """Calcula o vetor de deslocamentos da estrutura."""     
        
        self.u = np.linalg.inv(self.k01)@ self.fnos

    def calcula_solicitacoes_internas_nodais(self) -> None:
        """Calcula as forças no sistema local para as barras = solicitações."""
        
        for barra in self.barras:
            barra.monta_q()
            barra.calcula_r()

            # Monta o vetor u para as barras
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]

            #Calcula o vetor f para as barras

            barra.f = barra.k @ barra.u + barra.fep

            #Transforma o vetor f para o sistema local para as barras = solicitacoes
                 
            barra.fl = barra.r @ barra.f 

            barra.calcula_vetor_ASL()

    def calcula_reacoes(self) -> None:
        """Calcula um vetor com as reações da estrutura."""      
        
        for barra in self.barras:
            barra.monta_q()

            #Monta o vetor u para as barras
            
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]
            
            #Calcula o vetor f para as barras       
            
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
    """Representa as propriedades e parâmetros referente ao objeto carregamento distribuído."""
    
    def __init__(self, a:float, lw:float, w1:float, w2:float, barra:Barra) -> None:
        """
        Instancia um objeto da classe Carregamento_distribuido


        Args:

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
        
    def calcula_fepl(self) -> None:
        """Calcula o vetor de forças de engastamento perfeito no sistema local da barra."""

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

    def calcula_fep(self) -> None:
        """Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra."""        

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl

class Carregamento_pontual:
    """Representa as propriedades e parâmetros referente ao objeto carregamento pontual."""
   
    def __init__(self, a:float, Px:float, Py:float, barra:Barra) -> None:
        """
        Instancia um objeto da classe Carregamento_pontual


        Args:

            a: distância do nó inicial da barra até o ponto de aplicação da força
            Px: componente x da força
            Py: componente y da força
            barra: objeto barra
        """
        
        self.a = a
        self.Px = Px
        self.Py = Py
        self.barra = barra

    def calcula_fepl(self) -> None:
        """Calcula o vetor de forças de engastamento perfeito no sistema local da barra."""
        
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

    def calcula_fep(self) -> None:
        """Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra."""

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl

class Modelo:
    """Representa as propriedades e parâmetros referente ao objeto modelo."""

    def __init__(self, bw:float, h:float, l:float, lf:float, a:float, b:float, c:float, ap1:float, ap2:float, w:float, fck:float, fyk:float, fywk:float, modelo:int, α:float = 90.0, θ:float = 30.0, αE:float = 1.0, d:float = None, d_linha:float = None) -> None:
        """
        Instancia um objeto da classe Modelo


        Args:

            bw: largura viga
            h: altura da viga
            d: altura útil da armadura tracionada
            d_linha: altura útil da armadura comprimida 
            l: comprimento da viga
            lf: comprimento até centro da abertura
            a: altura da abertura
            b: largura da abertura
            c: espessura do banzo inferior na regiao da abertura
            ap1: largura pilar esquerdo
            ap2: largura pilar direito
            w: modulo carregamento distribuido
            fck: resistência característica do concreto [MPa]
            fyk: resistência característica do aço [MPa]
            fywk: resistência característica do aço para armadura transversal [MPa]
            αE: coeficiente de ajuste do módulo de elasticidade
            modelo: modelo de cálculo da armadura transversal (modelo 1 ou modelo 2)
            α: inclinação dos estribos (α deve estar entre 45° e 90°)
            θ: inclinação da biela comprimida de concreto (necessário para modelo 2, θ deve estar entre 30° e 45°)
        """

        self.bw = bw
        self.h = h
        self.d = d
        self.d_linha = d_linha
        self.l = l
        self.lf = lf
        self.a = a
        self.b = b
        self.c = c
        self.ap1 = ap1
        self.ap2 = ap2
        self.w = w
        self.fck = fck
        self.fyk = fyk
        self.fywk = fywk
        self.αE = αE
        self.modelo = modelo
        self.α = α
        self.θ = θ
        self.nos = None
        self.barras = None

        #Definição dos nós
        
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

        #Aplicação das restrições nodais
        
        no1.Tx = True
        no1.Ty = True
        no1.Rz = True
        no10.Tx = True
        no10.Ty = True
        no10.Rz = True

        # Definição dos materiais para uma seção de concreto armado
        concreto = Concreto(self.fck, self.αE)
        aco = Aco(self.fyk, self.fywk)

        # Definição da seção transversal

        secao1 = Secao("retangular", self.bw, self.h, concreto.Eci, self.h - 3.5, 3.5)
        secao2 = Secao("retangular", self.bw, self.c, concreto.Eci, self.c - 3.5, 3.5)
        secao3 = Secao("retangular", self.bw, self.h-self.a-self.c, concreto.Eci, self.h-self.a-self.c - 3.5, 3.5)
        secao4 = Secao("retangular", self.bw, self.h, concreto.Eci, self.h - 3.5, 3.5)

        secao4.I = secao4.I * 1000
        secao4.A = secao4.A * 1000

        #Definição das barras e propriedades
        #Ordem dos nós na definição das barras deve ser sempre da esquerda para direita ou de baixo para cima

        barra1 = Barra(1,no1,no2,secao1,concreto,aco)
        barra2 = Barra(2,no2,no3,secao1,concreto,aco)
        barra3 = Barra(3,no4,no7,secao2,concreto,aco)
        barra4 = Barra(4,no5,no8,secao3,concreto,aco)
        barra5 = Barra(5,no6,no9,secao1,concreto,aco)
        barra6 = Barra(6,no9,no10,secao1,concreto,aco)
        barra7 = Barra(7,no4,no3,secao4,concreto,aco)
        barra8 = Barra(8,no3,no5,secao4,concreto,aco)
        barra9 = Barra(9,no7,no6,secao4,concreto,aco)
        barra10 = Barra(10,no6,no8,secao4,concreto,aco) 
        self.barras = [barra1, barra2, barra3, barra4, barra5, barra6, barra7, barra8, barra9, barra10]

        #Definição dos carregamentos nas barras
        
        carregamento1 = Carregamento_distribuido(0.0,self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra2)
        carregamento2 = Carregamento_distribuido(0.0,self.b+self.bw,self.w,self.w,barra4)
        carregamento3 = Carregamento_distribuido(0.0,self.l-self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra5)
        
        #Calcula as forças de engastamento perfeito
    
        carregamento1.calcula_fepl()
        carregamento1.calcula_fep()
        carregamento2.calcula_fepl()
        carregamento2.calcula_fep()
        carregamento3.calcula_fepl()
        carregamento3.calcula_fep()
        
        #Define a estrutura, calcula os resultados de deslocamentos, solicitações internas e reações
    
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

        for barra in portico.barras:
            print("Asl barra" + str(barra.num))
            print_vetor(barra.ASL)

def print_matriz(m) -> None:
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


def print_vetor(v) -> None:
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












    





    




