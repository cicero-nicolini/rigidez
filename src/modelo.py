from .no import No
from .concreto import Concreto
from .aco import Aco
from .secao import Secao
from .barra import Barra
from .carregamento_distribuido import Carregamento_distribuido
from .estrutura import Estrutura
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from .visualizacao import print_vetor, print_matriz

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
        self.cord_x = [] 
        self.cord_y = []
        self.linhas = []
        self.x_asl = []
        self.y_asl = [] 
        self.linhas_asl = []

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
        
        no1.Tx = False
        no1.Ty = True
        no1.Rz = False
        no10.Tx = True
        no10.Ty = True
        no10.Rz = False

        # Definição dos materiais para uma seção de concreto armado
        concreto = Concreto(self.fck, self.αE)
        aco = Aco(self.fyk, self.fywk)

        # Definição da seção transversal

        d1 = self.h - 3.5
        d1_linha = 3.5

        d2 = self.c - 3.5
        d2_linha = 3.5

        d3 = self.h-self.a-self.c - 3.5
        d3_linha = 3.5

        d4 = self.h - 3.5
        d4_linha = 3.5

        secao1 = Secao("retangular", self.bw, self.h, concreto.Ecs, d1, d1_linha)
        secao2 = Secao("retangular", self.bw, self.c, concreto.Ecs, d2, d2_linha)
        secao3 = Secao("retangular", self.bw, self.h-self.a-self.c, concreto.Ecs, d3, d3_linha)
        secao4 = Secao("retangular", 100.0, 100.0, 100.0*concreto.Ecs, 0.0, 0.0)

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

        carregamento1 = Carregamento_distribuido(0.0,ap1/2,self.w,self.w,barra1)
        carregamento2 = Carregamento_distribuido(0.0,self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra2)
        carregamento3 = Carregamento_distribuido(0.0,self.b+self.bw,self.w,self.w,barra4)
        carregamento4 = Carregamento_distribuido(0.0,self.l-self.lf-(self.b/2)-(self.bw/2),self.w,self.w,barra5)
        carregamento5 = Carregamento_distribuido(0.0,ap2/2,self.w,self.w,barra6)

        #Calcula as forças de engastamento perfeito
    
        carregamento1.calcula_fepl()
        carregamento1.calcula_fep()
        carregamento2.calcula_fepl()
        carregamento2.calcula_fep()
        carregamento3.calcula_fepl()
        carregamento3.calcula_fep()
        carregamento4.calcula_fepl()
        carregamento4.calcula_fep()
        carregamento5.calcula_fepl()
        carregamento5.calcula_fep()
        
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

        portico.calcula_ASL_barra()
        
        for barra in portico.barras:
            
            print("Asl barra" + str(barra.num))
            print_vetor(barra.ASL)

        for barra in self.barras:
                    
            print("comprimento barra" + str(barra.num))
            print(barra.L) 

        for no in self.nos:
            
            print("no" + str(no.num))
            print(no.x)
            print(no.y)

    def desenha_area_armaduras(self) -> None:
        """Desenha gráfico com resultados de áreas de armaduras por nó da estrutura."""
        
        for no in self.nos:
            self.cord_y.append(no.y)  
            self.cord_x.append(no.x)

        for barra in self.barras:

            no_inicial = (barra.noi.x, barra.noi.y)    
            no_final = (barra.noj.x, barra.noj.y)

            self.linhas.append([no_inicial, no_final])

        #lista com barras onde serão mostrados os resultados para área de armaduras
        barras_asl = [self.barras[1], self.barras[2], self.barras[3], self.barras[4]]
              
        #nó 1 estrutura - barra 1
        self.x_asl.append(self.barras[0].noi.x)
        self.y_asl.append(self.barras[0].noi.y - self.barras[0].ASL[0])

        #nó 10 - barra 6
        self.x_asl.append(self.barras[5].noj.x)
        self.y_asl.append(self.barras[5].noj.y - self.barras[5].ASL[1])

        #linha do diagrama que liga valores asl referentes ao nó 1 e nó 2 da estrutura
        no1_asl = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASL[0])    
        no2_asl = (self.barras[0].noj.x, self.barras[0].noj.y - self.barras[0].ASL[1])
        self.linhas_asl.append([no1_asl, no2_asl])

        #linha do diagrama que liga valores asl referentes ao nó 9 e nó 10 da estrutura
        no9_asl = (self.barras[5].noi.x, self.barras[5].noi.y - self.barras[5].ASL[0])    
        no10_asl = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASL[1])
        self.linhas_asl.append([no9_asl, no10_asl])

        #pontos para formar a primeira linha vertical que delimita diagrama
        ponto1_origem_diagrama = (self.barras[0].noi.x, self.barras[0].noi.y)    
        ponto2_origem_diagrama = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASL[0])
        self.linhas_asl.append([ponto1_origem_diagrama, ponto2_origem_diagrama])

        #pontos para formar a ultima linha vertical que delimita diagrama
        ponto1_fim_diagrama = (self.barras[5].noj.x, self.barras[5].noj.y)    
        ponto2_fim_diagrama = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASL[1])
        self.linhas_asl.append([ponto1_fim_diagrama, ponto2_fim_diagrama])
        
        for barra in barras_asl:

            self.x_asl.append(barra.noi.x)
            self.y_asl.append(barra.noi.y - barra.ASL[0])
           
            self.x_asl.append(barra.noj.x)
            self.y_asl.append(barra.noj.y - barra.ASL[1])

            no1_linha_asl = (barra.noi.x, barra.noi.y - barra.ASL[0])    
            no2_linha_asl = (barra.noj.x, barra.noj.y - barra.ASL[1])
            self.linhas_asl.append([no1_linha_asl, no2_linha_asl])

            #pontos para formar linha vertical que delimita diagrama para cada barra
            
            ponto1_barra_origem_diagrama = (barra.noi.x, barra.noi.y) 
            ponto2_barra_origem_diagrama = (barra.noi.x, barra.noi.y - barra.ASL[0])
            self.linhas_asl.append([ponto1_barra_origem_diagrama, ponto2_barra_origem_diagrama])

            ponto1_barra_fim_diagrama = (barra.noj.x, barra.noj.y) 
            ponto2_barra_fim_diagrama = (barra.noj.x, barra.noj.y - barra.ASL[1])
            self.linhas_asl.append([ponto1_barra_fim_diagrama, ponto2_barra_fim_diagrama])
            
        lc = LineCollection(self.linhas, colors='dimgray', linewidths=1.5)
        lc2 = LineCollection(self.linhas_asl, colors='royalblue', linewidths=1)
        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        ax.add_collection(lc)
        ax.add_collection(lc2)
        ax.plot(self.cord_x, self.cord_y, 'ko', markersize=3)
        ax.plot(self.x_asl, self.y_asl, 'ro', markersize=3)
        ax.set_aspect(5)
        plt.show()  