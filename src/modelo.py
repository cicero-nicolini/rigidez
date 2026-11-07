from .cargas import normalizar_cargas, analisar_combinacoes, selecionar_esforcos
import math
from .vinculacoes import aplicar_vinculacao, validar_vinculacao
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

    def __init__(self, bw:float, h:float, l:float, lf:float, a:float, b:float, c:float, ap1:float, ap2:float, w:float = None, fck:float = None, fyk:float = None, fywk:float = None, modelo:int = 1, α:float = 90.0, θ:float = 30.0, αE:float = 1.0, cob:float = 3.5, espacamento_maximo:float = None, vinculacao_esquerda:int = 1, vinculacao_direita:int = 2, cargas=None) -> None:
        """
        Instancia um objeto da classe Modelo


        Args:

            bw: largura viga
            h: altura da viga
            cob: distância da face do concreto ao centroide da armadura longitudinal (cm)
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
            vinculacao_esquerda: classe de apoio esquerdo (1, 2 ou 3).
            vinculacao_direita: classe de apoio direito (1, 2 ou 3).
            espacamento_maximo: comprimento máximo de cada elemento, em cm.
                Se None, mantém oito elementos por trecho horizontal principal
                e um elemento nos demais trechos (malha original).
                As barras de seção fictícia nunca são subdivididas.
        """

        if any(v is None for v in (fck, fyk, fywk)):
            raise ValueError("Informe fck, fyk e fywk.")
        self.vinculacao_esquerda = validar_vinculacao(vinculacao_esquerda)
        self.vinculacao_direita = validar_vinculacao(vinculacao_direita)
        self.bw = bw
        self.h = h
        self.a = a
        self.b = b
        self.c = c
        self.cob = float(cob)
        if not math.isfinite(self.cob) or self.cob <= 0:
            raise ValueError("cob deve ser positivo e finito.")
        if 2 * self.cob >= min(self.h, self.c, self.h - self.a - self.c):
            raise ValueError("cob incompatível com a altura de alguma seção.")
        self.l = l
        self.lf = lf
        self.a = a
        self.b = b
        self.c = c
        self.ap1 = ap1
        self.ap2 = ap2
        self.cargas = normalizar_cargas(cargas, w)
        self.wk = sum(c["w"] for c in self.cargas)
        self.wd = sum(c["w"] * c["gamaf"] for c in self.cargas)
        self.w = self.wk  # intensidade característica para visualização
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

        # Pontos obrigatórios: apoios, encontros dos banzos e montantes.
        x_esq = self.ap1 / 2.0 + self.lf - (self.b + self.bw) / 2.0
        x_dir = self.ap1 / 2.0 + self.lf + (self.b + self.bw) / 2.0
        y_inf = -(self.h - self.c) / 2.0
        y_sup = (self.a + self.c) / 2.0
        coordenadas = [
            (0.0, 0.0), (self.ap1 / 2.0, 0.0),
            (x_esq, 0.0), (x_esq, y_inf), (x_esq, y_sup),
            (x_dir, 0.0), (x_dir, y_inf), (x_dir, y_sup),
            (self.ap1 / 2.0 + self.l, 0.0),
            (self.ap1 / 2.0 + self.l + self.ap2 / 2.0, 0.0),
        ]
        if espacamento_maximo is not None:
            if not math.isfinite(espacamento_maximo) or espacamento_maximo <= 0:
                raise ValueError('espacamento_maximo deve ser positivo e finito.')
        self.espacamento_maximo = espacamento_maximo
        self.nos = [No(i, float(x), float(y))
                    for i, (x, y) in enumerate(coordenadas, start=1)]
        aplicar_vinculacao(self.nos[0], self.vinculacao_esquerda)
        aplicar_vinculacao(self.nos[9], self.vinculacao_direita)

        # Definição dos materiais para uma seção de concreto armado
        concreto = Concreto(self.fck, self.αE)
        aco = Aco(self.fyk, self.fywk)

        # Definição da seção transversal

        d1 = self.h - self.cob
        d1_linha = self.cob

        d2 = self.c - self.cob
        d2_linha = self.cob

        d3 = self.h-self.a-self.c - self.cob
        d3_linha = self.cob


        secao1 = Secao("retangular", self.bw, self.h, concreto.Ecs, d1, d1_linha)
        secao2 = Secao("retangular", self.bw, self.c, concreto.Ecs, d2, d2_linha)
        secao3 = Secao("retangular", self.bw, self.h-self.a-self.c, concreto.Ecs, d3, d3_linha)
        secao4 = Secao("retangular", 100.0, 100.0, 100.0*concreto.Ecs, 0.0, 0.0)

        # Trechos de discretização. Cada trecho recebe elementos de
        # comprimento idêntico, sempre <= espacamento_maximo.
        # Mantém a orientação das barras: esquerda-direita ou baixo-cima.
        trechos = [
            (1, 2, secao1, True, 1),
            (2, 3, secao1, True, 8),
            (4, 7, secao2, False, 8),
            (5, 8, secao3, True, 8),
            (6, 9, secao1, True, 8),
            (9, 10, secao1, True, 1),
            (4, 3, secao4, False, 1),
            (3, 5, secao4, False, 1),
            (7, 6, secao4, False, 1),
            (6, 8, secao4, False, 1),
        ]
        self.barras = []
        self.carregamentos = []
        for indice_trecho, (ni, nj, secao, carregado, n_original) in enumerate(trechos):
            no_inicio = self.nos[ni - 1]
            no_fim = self.nos[nj - 1]
            comprimento = math.hypot(no_fim.x - no_inicio.x,
                                     no_fim.y - no_inicio.y)
            if comprimento <= 0:
                raise ValueError(f'Trecho {ni}-{nj} com comprimento nulo.')
            # Os montantes com seção fictícia (secao4) são elementos
            # de ligação rígida e não devem ser discretizados.
            # O espaçamento máximo aplica-se somente às barras reais.
            if secao is secao4:
                n_elementos = 1
            elif espacamento_maximo is None:
                n_elementos = n_original
            else:
                n_elementos = max(1, math.ceil(comprimento / espacamento_maximo))
            nos_trecho = [no_inicio]
            for k in range(1, n_elementos):
                t = k / n_elementos
                no = No(len(self.nos) + 1,
                        no_inicio.x + t * (no_fim.x - no_inicio.x),
                        no_inicio.y + t * (no_fim.y - no_inicio.y))
                self.nos.append(no)
                nos_trecho.append(no)
            nos_trecho.append(no_fim)
            for no_i, no_j in zip(nos_trecho[:-1], nos_trecho[1:]):
                barra = Barra(len(self.barras) + 1, no_i, no_j,
                              secao, concreto, aco)
                barra.ficticia = (secao is secao4)
                self.barras.append(barra)
                barra.carregada = bool(carregado)
                if carregado:
                    carga = Carregamento_distribuido(
                        0.0, barra.L, self.w, self.w, barra)
                    self.carregamentos.append(carga)
                    carga.calcula_fepl()
                    carga.calcula_fep()

        self.portico = Estrutura(self.nos, self.barras)
        self.portico.monta_k()
        self.portico.monta_k01()
        analisar_combinacoes(self)

    def selecionar_esforcos(self, tipo="calculo"):
        """Seleciona os resultados característicos ou majorados."""
        selecionar_esforcos(self, tipo)

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