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
        self.x_asw_s = []
        self.y_asw_s = [] 
        self.linhas_asw_s = []

        #Definição dos nós
        
        no1 = No(1, 0.0, 0.0)
        no2 = No(2, self.ap1/2.0, 0.0)
        no11 = No(11, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(1.0/8.0)), 0.0)
        no12 = No(12, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(2.0/8.0)), 0.0)
        no13 = No(13, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(3.0/8.0)), 0.0)
        no14 = No(14, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(4.0/8.0)), 0.0)
        no15 = No(15, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(5.0/8.0)), 0.0)
        no16 = No(16, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(6.0/8.0)), 0.0)
        no17 = No(17, self.ap1/2.0 + ((self.lf - self.b/2.0 - self.bw/2.0)*(7.0/8.0)), 0.0)
        no3 = No(3, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, 0.0)
        no4 = No(4, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, -(self.h-self.c)/2.0)
        no18 = No(18, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(1.0/8.0)), -(self.h-self.c)/2.0)
        no19 = No(19, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(2.0/8.0)), -(self.h-self.c)/2.0)
        no20 = No(20, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(3.0/8.0)), -(self.h-self.c)/2.0)
        no21 = No(21, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(4.0/8.0)), -(self.h-self.c)/2.0)
        no22 = No(22, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(5.0/8.0)), -(self.h-self.c)/2.0)
        no23 = No(23, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(6.0/8.0)), -(self.h-self.c)/2.0)
        no24 = No(24, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(7.0/8.0)), -(self.h-self.c)/2.0)
        no7 = No(7, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, -(self.h-self.c)/2.0)
        no5 = No(5, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0, (self.a+self.c)/2.0)
        no25 = No(25, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(1.0/8.0)), (self.a+self.c)/2.0)
        no26 = No(26, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(2.0/8.0)), (self.a+self.c)/2.0)
        no27 = No(27, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(3.0/8.0)), (self.a+self.c)/2.0)
        no28 = No(28, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(4.0/8.0)), (self.a+self.c)/2.0)
        no29 = No(29, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(5.0/8.0)), (self.a+self.c)/2.0)
        no30 = No(30, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(6.0/8.0)), (self.a+self.c)/2.0)
        no31 = No(31, self.ap1/2.0 + self.lf - self.b/2.0 - self.bw/2.0 + ((self.b + self.bw)*(7.0/8.0)), (self.a+self.c)/2.0)
        no8 = No(8, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, (self.a+self.c)/2.0)
        no6 = No(6, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0, 0.0)
        no32 = No(32, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(1.0/8.0)), 0.0)
        no33 = No(33, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(2.0/8.0)), 0.0)
        no34 = No(34, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(3.0/8.0)), 0.0)
        no35 = No(35, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(4.0/8.0)), 0.0)
        no36 = No(36, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(5.0/8.0)), 0.0)
        no37 = No(37, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(6.0/8.0)), 0.0)
        no38 = No(38, self.ap1/2.0 + self.lf + self.b/2.0 + self.bw/2.0 + ((self.l - self.lf - self.b/2.0 - self.bw/2.0)*(7.0/8.0)), 0.0)
        no9 = No(9, self.ap1/2.0 + self.l, 0.0)
        no10 = No(10, self.ap1/2.0 + self.l + self.ap2/2.0, 0.0)
        self.nos = [no1, no2, no3, no4, no5, no6, no7, no8, no9, no10, no11, no12, no13,
                    no14, no15, no16, no17, no18, no19, no20, no21, no22, no23, no24, no25,
                    no26, no27, no28, no29, no30, no31, no32, no33, no34, no35, no36, no37, no38]

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
        barra2 = Barra(2,no2,no11,secao1,concreto,aco)
        barra11 = Barra(11,no11,no12,secao1,concreto,aco)
        barra12 = Barra(12,no12,no13,secao1,concreto,aco)
        barra13 = Barra(13,no13,no14,secao1,concreto,aco)
        barra14 = Barra(14,no14,no15,secao1,concreto,aco)
        barra15 = Barra(15,no15,no16,secao1,concreto,aco)
        barra16 = Barra(16,no16,no17,secao1,concreto,aco)
        barra17 = Barra(17,no17,no3,secao1,concreto,aco)
        barra3 = Barra(3,no4,no18,secao2,concreto,aco)
        barra18 = Barra(18,no18,no19,secao2,concreto,aco)
        barra19 = Barra(19,no19,no20,secao2,concreto,aco)
        barra20 = Barra(20,no20,no21,secao2,concreto,aco)
        barra21 = Barra(21,no21,no22,secao2,concreto,aco)
        barra22 = Barra(22,no22,no23,secao2,concreto,aco)
        barra23 = Barra(23,no23,no24,secao2,concreto,aco)
        barra24 = Barra(24,no24,no7,secao2,concreto,aco)
        barra4 = Barra(4,no5,no25,secao3,concreto,aco)
        barra25 = Barra(25,no25,no26,secao3,concreto,aco)
        barra26 = Barra(26,no26,no27,secao3,concreto,aco)
        barra27 = Barra(27,no27,no28,secao3,concreto,aco)
        barra28 = Barra(28,no28,no29,secao3,concreto,aco)
        barra29 = Barra(29,no29,no30,secao3,concreto,aco)
        barra30 = Barra(30,no30,no31,secao3,concreto,aco)
        barra31 = Barra(31,no31,no8,secao3,concreto,aco)
        barra5 = Barra(5,no6,no32,secao1,concreto,aco)
        barra32 = Barra(32,no32,no33,secao1,concreto,aco)
        barra33 = Barra(33,no33,no34,secao1,concreto,aco)
        barra34 = Barra(34,no34,no35,secao1,concreto,aco)
        barra35 = Barra(35,no35,no36,secao1,concreto,aco)
        barra36 = Barra(36,no36,no37,secao1,concreto,aco)
        barra37 = Barra(37,no37,no38,secao1,concreto,aco)
        barra38 = Barra(38,no38,no9,secao1,concreto,aco)
        barra6 = Barra(6,no9,no10,secao1,concreto,aco)
        barra7 = Barra(7,no4,no3,secao4,concreto,aco)
        barra8 = Barra(8,no3,no5,secao4,concreto,aco)
        barra9 = Barra(9,no7,no6,secao4,concreto,aco)
        barra10 = Barra(10,no6,no8,secao4,concreto,aco) 
        self.barras = [barra1, barra2, barra3, barra4, barra5, barra6, barra7, barra8, barra9, 
                       barra10, barra11, barra12, barra13, barra14, barra15, barra16, barra17,
                       barra18, barra19, barra20, barra21, barra22, barra23, barra24, barra25,
                       barra26, barra27, barra28, barra29, barra30, barra31, barra32, barra33,
                       barra34, barra35, barra36, barra37, barra38]

        #Definição dos carregamentos nas barras

        carregamento1 = Carregamento_distribuido(0.0,ap1/2,self.w,self.w,barra1)
        carregamento2 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra2)
        carregamento6 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra11)
        carregamento7 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra12)
        carregamento8 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra13)
        carregamento9 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra14)
        carregamento10 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra15)
        carregamento11 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra16)
        carregamento12 = Carregamento_distribuido(0.0,((self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra17)
        carregamento3 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra4)
        carregamento13 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra25)
        carregamento14 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra26)
        carregamento15 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra27)
        carregamento16 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra28)
        carregamento17 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra29)
        carregamento18 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra30)
        carregamento19 = Carregamento_distribuido(0.0,((self.b+self.bw)/8.0),self.w,self.w,barra31)
        carregamento4 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra5)
        carregamento20 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra32)
        carregamento21 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra33)
        carregamento22 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra34)
        carregamento23 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra35)
        carregamento24 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra36)
        carregamento25 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra37)
        carregamento26 = Carregamento_distribuido(0.0,((self.l-self.lf-(self.b/2)-(self.bw/2))/8.0),self.w,self.w,barra38)
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
        carregamento6.calcula_fepl()
        carregamento6.calcula_fep()
        carregamento7.calcula_fepl()
        carregamento7.calcula_fep()
        carregamento8.calcula_fepl()
        carregamento8.calcula_fep()
        carregamento9.calcula_fepl()
        carregamento9.calcula_fep()
        carregamento10.calcula_fepl()
        carregamento10.calcula_fep()
        carregamento11.calcula_fepl()
        carregamento11.calcula_fep()
        carregamento12.calcula_fepl()
        carregamento12.calcula_fep()
        carregamento13.calcula_fepl()
        carregamento13.calcula_fep()
        carregamento14.calcula_fepl()
        carregamento14.calcula_fep()
        carregamento15.calcula_fepl()
        carregamento15.calcula_fep()
        carregamento16.calcula_fepl()
        carregamento16.calcula_fep()
        carregamento17.calcula_fepl()
        carregamento17.calcula_fep()
        carregamento18.calcula_fepl()
        carregamento18.calcula_fep()
        carregamento19.calcula_fepl()
        carregamento19.calcula_fep()
        carregamento20.calcula_fepl()
        carregamento20.calcula_fep()
        carregamento21.calcula_fepl()
        carregamento21.calcula_fep()
        carregamento22.calcula_fepl()
        carregamento22.calcula_fep()
        carregamento23.calcula_fepl()
        carregamento23.calcula_fep()
        carregamento24.calcula_fepl()
        carregamento24.calcula_fep()
        carregamento25.calcula_fepl()
        carregamento25.calcula_fep()
        carregamento26.calcula_fepl()
        carregamento26.calcula_fep()
        
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
        print("solicitacoes barra11")
        print_vetor(portico.barras[10].fl)
        print("solicitacoes barra12")
        print_vetor(portico.barras[11].fl)
        print("solicitacoes barra13")
        print_vetor(portico.barras[12].fl)
        print("solicitacoes barra14")
        print_vetor(portico.barras[13].fl)
        print("solicitacoes barra15")
        print_vetor(portico.barras[14].fl)
        print("solicitacoes barra16")
        print_vetor(portico.barras[15].fl)
        print("solicitacoes barra17")
        print_vetor(portico.barras[16].fl)
        print("solicitacoes barra18")
        print_vetor(portico.barras[17].fl)
        print("solicitacoes barra19")
        print_vetor(portico.barras[18].fl)
        print("solicitacoes barra20")
        print_vetor(portico.barras[19].fl)
        print("solicitacoes barra21")
        print_vetor(portico.barras[20].fl)
        print("solicitacoes barra22")
        print_vetor(portico.barras[21].fl)
        print("solicitacoes barra23")
        print_vetor(portico.barras[22].fl)
        print("solicitacoes barra24")
        print_vetor(portico.barras[23].fl)
        print("solicitacoes barra25")
        print_vetor(portico.barras[24].fl)
        print("solicitacoes barra26")
        print_vetor(portico.barras[25].fl)
        print("solicitacoes barra27")
        print_vetor(portico.barras[26].fl)
        print("solicitacoes barra28")
        print_vetor(portico.barras[27].fl)
        print("solicitacoes barra29")
        print_vetor(portico.barras[28].fl)
        print("solicitacoes barra30")
        print_vetor(portico.barras[29].fl)
        print("solicitacoes barra31")
        print_vetor(portico.barras[30].fl)
        print("solicitacoes barra32")
        print_vetor(portico.barras[31].fl)
        print("solicitacoes barra33")
        print_vetor(portico.barras[32].fl)
        print("solicitacoes barra34")
        print_vetor(portico.barras[33].fl)
        print("solicitacoes barra35")
        print_vetor(portico.barras[34].fl)
        print("solicitacoes barra36")
        print_vetor(portico.barras[35].fl)
        print("solicitacoes barra37")
        print_vetor(portico.barras[36].fl)
        print("solicitacoes barra38")
        print_vetor(portico.barras[37].fl)

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
        barras_as = [self.barras[1],self.barras[2], self.barras[3], self.barras[4],self.barras[10],
                      self.barras[11],self.barras[12], self.barras[13], self.barras[14],self.barras[15],
                      self.barras[16],self.barras[17], self.barras[18], self.barras[19],self.barras[20],
                      self.barras[21],self.barras[22], self.barras[23], self.barras[24],self.barras[25],
                      self.barras[26],self.barras[27], self.barras[28], self.barras[29],self.barras[30],
                      self.barras[31],self.barras[32], self.barras[33], self.barras[34],self.barras[35],
                      self.barras[36],self.barras[37]]
              
        #nó 1 estrutura - barra 1
        self.x_asl.append(self.barras[0].noi.x)
        self.y_asl.append(self.barras[0].noi.y - self.barras[0].ASL[0])

        self.x_asw_s.append(self.barras[0].noi.x)
        self.y_asw_s.append(self.barras[0].noi.y - self.barras[0].ASW_S[0])

        #nó 10 - barra 6
        self.x_asl.append(self.barras[5].noj.x)
        self.y_asl.append(self.barras[5].noj.y - self.barras[5].ASL[1])

        self.x_asw_s.append(self.barras[5].noj.x)
        self.y_asw_s.append(self.barras[5].noj.y - self.barras[5].ASW_S[1])

        #linha do diagrama que liga valores asl referentes ao nó 1 e nó 2 da estrutura
        no1_asl = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASL[0])    
        no2_asl = (self.barras[0].noj.x, self.barras[0].noj.y - self.barras[0].ASL[1])
        self.linhas_asl.append([no1_asl, no2_asl])

        no1_asw_s = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASW_S[0])    
        no2_asw_s = (self.barras[0].noj.x, self.barras[0].noj.y - self.barras[0].ASW_S[1])
        self.linhas_asl.append([no1_asw_s, no2_asw_s])

        #linha do diagrama que liga valores asl referentes ao nó 9 e nó 10 da estrutura
        no9_asl = (self.barras[5].noi.x, self.barras[5].noi.y - self.barras[5].ASL[0])    
        no10_asl = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASL[1])
        self.linhas_asl.append([no9_asl, no10_asl])

        no9_asw_s = (self.barras[5].noi.x, self.barras[5].noi.y - self.barras[5].ASW_S[0])    
        no10_asw_s = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASW_S[1])
        self.linhas_asl.append([no9_asw_s, no10_asw_s])

        #pontos para formar a primeira linha vertical que delimita diagrama
        ponto1_origem_diagrama_asl = (self.barras[0].noi.x, self.barras[0].noi.y)    
        ponto2_origem_diagrama_asl = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASL[0])
        self.linhas_asl.append([ponto1_origem_diagrama_asl, ponto2_origem_diagrama_asl])

        ponto1_origem_diagrama_asw_s = (self.barras[0].noi.x, self.barras[0].noi.y)    
        ponto2_origem_diagrama_asw_s = (self.barras[0].noi.x, self.barras[0].noi.y - self.barras[0].ASL[1])
        self.linhas_asl.append([ponto1_origem_diagrama_asw_s, ponto2_origem_diagrama_asw_s])

        #pontos para formar a ultima linha vertical que delimita diagrama
        ponto1_fim_diagrama_asl = (self.barras[5].noj.x, self.barras[5].noj.y)    
        ponto2_fim_diagrama_asl = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASL[1])
        self.linhas_asl.append([ponto1_fim_diagrama_asl, ponto2_fim_diagrama_asl])

        ponto1_fim_diagrama_asw_s = (self.barras[5].noj.x, self.barras[5].noj.y)    
        ponto2_fim_diagrama_asw_s = (self.barras[5].noj.x, self.barras[5].noj.y - self.barras[5].ASW_S[1])
        self.linhas_asl.append([ponto1_fim_diagrama_asw_s, ponto2_fim_diagrama_asw_s])
        
        for barra in barras_as:

            self.x_asl.append(barra.noi.x)
            self.y_asl.append(barra.noi.y - barra.ASL[0])

            self.x_asw_s.append(barra.noi.x)
            self.y_asw_s.append(barra.noi.y - barra.ASW_S[0])
           
            self.x_asl.append(barra.noj.x)
            self.y_asl.append(barra.noj.y - barra.ASL[1])

            self.x_asw_s.append(barra.noj.x)
            self.y_asw_s.append(barra.noj.y - barra.ASW_S[1])

            no1_linha_asl = (barra.noi.x, barra.noi.y - barra.ASL[0])    
            no2_linha_asl = (barra.noj.x, barra.noj.y - barra.ASL[1])
            self.linhas_asl.append([no1_linha_asl, no2_linha_asl])

            no1_linha_asw_s = (barra.noi.x, barra.noi.y - barra.ASW_S[0])    
            no2_linha_asw_s = (barra.noj.x, barra.noj.y - barra.ASW_S[1])
            self.linhas_asw_s.append([no1_linha_asw_s, no2_linha_asw_s])

            #pontos para formar linha vertical que delimita diagrama para cada barra
            
            ponto1_barra_origem_diagrama_asl = (barra.noi.x, barra.noi.y) 
            ponto2_barra_origem_diagrama_asl = (barra.noi.x, barra.noi.y - barra.ASL[0])
            self.linhas_asl.append([ponto1_barra_origem_diagrama_asl, ponto2_barra_origem_diagrama_asl])

            ponto1_barra_origem_diagrama_asw_s = (barra.noi.x, barra.noi.y) 
            ponto2_barra_origem_diagrama_asw_s = (barra.noi.x, barra.noi.y - barra.ASW_S[0])
            self.linhas_asw_s.append([ponto1_barra_origem_diagrama_asw_s, ponto2_barra_origem_diagrama_asw_s])

            ponto1_barra_fim_diagrama_asl = (barra.noj.x, barra.noj.y) 
            ponto2_barra_fim_diagrama_asl = (barra.noj.x, barra.noj.y - barra.ASL[1])
            self.linhas_asl.append([ponto1_barra_fim_diagrama_asl, ponto2_barra_fim_diagrama_asl])

            ponto1_barra_fim_diagrama_asw_s = (barra.noj.x, barra.noj.y) 
            ponto2_barra_fim_diagrama_asw_s = (barra.noj.x, barra.noj.y - barra.ASW_S[1])
            self.linhas_asw_s.append([ponto1_barra_fim_diagrama_asw_s, ponto2_barra_fim_diagrama_asw_s])
            
        lc = LineCollection(
            self.linhas,
            colors='dimgray',
            linewidths=1.5
        )

        lc2 = LineCollection(
            self.linhas_asl,
            colors='#2563A6',
            linewidths=1.2
        )

        lc3 = LineCollection(
            self.linhas_asw_s,
            colors='#D97706',
            linewidths=1.2
        )

        fig, ax = plt.subplots()

        ax.add_collection(lc)
        ax.add_collection(lc2)
        ax.add_collection(lc3)

        ax.plot(self.cord_x, self.cord_y, 'ko', markersize=3)
        ax.plot(self.x_asl, self.y_asl, 'o', color='#2563A6', markersize=3)
        ax.plot(self.x_asw_s, self.y_asw_s, 'o', color='#D97706', markersize=3)

        ax.set_aspect(4)
        plt.show() 