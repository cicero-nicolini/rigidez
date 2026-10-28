
import numpy as np
import math
from .secao import Secao
from .no import No
from .concreto import Concreto
from .aco import Aco

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

    def calcula_vetor_ASW_S(self) -> None:
        """
        Calcula vetor ASW_S que resulta na área de armadura transversal em cm²/m.
        
        Raises
            ValueError: Se a força cortante solicitante de cálculo (Vsd) for superior
            à resistência das bielas comprimidas de concreto (Vrd2).
        """

        α = math.radians(self.α)
        θ = math.radians(self.θ)

        # Definindo o coeficiente majorador de esforços
        γf = 1.4

        # Guardando as solicitações (M, N) dos nós extremos da barra
        M = np.array([self.fl[2], self.fl[5]])
        N = np.array([self.fl[0], self.fl[3]])
        V = np.array([self.fl[1], self.fl[4]])

        Vsd = abs(V*γf)

        for i in range(len(V)):

            if self.modelo == 1:
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

            if self.modelo == 2:
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

            if N[0] < 0 and N[1] > 0:
            #if n > 0:
                # Flexo-tração
                e0 = abs(M/N) # calculo da excentricidade da carga
                
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

            if N[0] > 0 and N[1] < 0:
            #if n < 0:
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
       