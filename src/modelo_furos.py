"""Modelo de viga com múltiplas aberturas retangulares.

Mantém a formulação de pórtico plano, os banzos e as ligações rígidas
fictícias do modelo de abertura única.
"""
import math
from .furo import Furo
from .no import No
from .concreto import Concreto
from .aco import Aco
from .secao import Secao
from .barra import Barra
from .carregamento_distribuido import Carregamento_distribuido
from .estrutura import Estrutura
from .visualizacao import print_vetor, print_matriz


class ModeloFuros:
    """Viga com zero, uma ou várias aberturas; dimensões em cm."""

    def __init__(self, bw, h, l, furos, ap1, ap2, w, fck, fyk, fywk,
                 modelo, α=90.0, θ=30.0, αE=1.0, d=None, d_linha=None,
                 espacamento_maximo=20.0, imprimir_resultados=False):
        self.bw, self.h, self.l = float(bw), float(h), float(l)
        self.ap1, self.ap2, self.w = float(ap1), float(ap2), float(w)
        self.fck, self.fyk, self.fywk = fck, fyk, fywk
        self.modelo, self.α, self.θ, self.αE = modelo, α, θ, αE
        self.d, self.d_linha = d, d_linha
        self.espacamento_maximo = espacamento_maximo
        self.furos = sorted(list(furos) if furos is not None else [],
                            key=lambda f: f.lf)
        if any(not isinstance(f, Furo) for f in self.furos):
            raise TypeError('Todos os elementos de furos devem ser instâncias de Furo.')
        if min(self.bw, self.h, self.l, self.ap1, self.ap2) <= 0:
            raise ValueError('As dimensões da viga e dos pilares devem ser positivas.')
        if espacamento_maximo is not None and (
                not math.isfinite(espacamento_maximo) or espacamento_maximo <= 0):
            raise ValueError('espacamento_maximo deve ser positivo e finito.')
        for f in self.furos:
            f.validar_na_viga(self.l, self.h, self.bw)
        for anterior, seguinte in zip(self.furos, self.furos[1:]):
            if anterior.sobrepoe(seguinte, folga=self.bw):
                raise ValueError('As zonas de transição dos furos se sobrepõem.')

        self.nos, self.barras, self.carregamentos = [], [], []
        self.cord_x, self.cord_y, self.linhas = [], [], []
        self.x_asl, self.y_asl, self.linhas_asl = [], [], []
        self.x_asw_s, self.y_asw_s, self.linhas_asw_s = [], [], []

        concreto = Concreto(self.fck, self.αE)
        aco = Aco(self.fyk, self.fywk)
        def secao(altura):
            return Secao('retangular', self.bw, altura, concreto.Ecs,
                         altura - 3.5, 3.5)
        secao_cheia = secao(self.h)
        secao_rigida = Secao('retangular', 100.0, 100.0,
                             100.0 * concreto.Ecs, 0.0, 0.0)

        def no(x, y):
            n = No(len(self.nos) + 1, float(x), float(y))
            self.nos.append(n)
            return n

        def trecho(ni, nj, sec, carregado=False, ficticia=False):
            L = math.hypot(nj.x - ni.x, nj.y - ni.y)
            if L <= 1e-10:
                raise ValueError('Trecho de comprimento nulo.')
            quantidade = (1 if ficticia else
                          (8 if espacamento_maximo is None else
                           max(1, math.ceil(L / espacamento_maximo))))
            nos_trecho = [ni]
            for k in range(1, quantidade):
                t = k / quantidade
                nos_trecho.append(no(ni.x + t * (nj.x - ni.x),
                                     ni.y + t * (nj.y - ni.y)))
            nos_trecho.append(nj)
            for inicio, fim in zip(nos_trecho, nos_trecho[1:]):
                barra = Barra(len(self.barras) + 1, inicio, fim, sec, concreto, aco)
                barra.ficticia = ficticia
                self.barras.append(barra)
                if carregado:
                    carga = Carregamento_distribuido(0.0, barra.L,
                                                     self.w, self.w, barra)
                    self.carregamentos.append(carga)
                    carga.calcula_fepl()
                    carga.calcula_fep()

        apoio_esq = no(0.0, 0.0)
        apoio_esq.Ty = True
        inicio_viga = no(self.ap1 / 2.0, 0.0)
        trecho(apoio_esq, inicio_viga, secao_cheia, carregado=True)
        atual = inicio_viga
        for f in self.furos:
            xesq = self.ap1 / 2.0 + f.inicio - self.bw / 2.0
            xdir = self.ap1 / 2.0 + f.fim + self.bw / 2.0
            # Eixos dos banzos: mesma convenção do modelo original.
            yinf = -(self.h - f.c) / 2.0
            ysup = (f.a + f.c) / 2.0
            esquerda = no(xesq, 0.0)
            inferior_esq = no(xesq, yinf)
            superior_esq = no(xesq, ysup)
            direita = no(xdir, 0.0)
            inferior_dir = no(xdir, yinf)
            superior_dir = no(xdir, ysup)
            trecho(atual, esquerda, secao_cheia, carregado=True)
            trecho(inferior_esq, inferior_dir, secao(f.c))
            trecho(superior_esq, superior_dir, secao(self.h - f.a - f.c),
                   carregado=True)
            # Ligações rígidas: nunca discretizadas e sem dimensionamento.
            trecho(inferior_esq, esquerda, secao_rigida, ficticia=True)
            trecho(esquerda, superior_esq, secao_rigida, ficticia=True)
            trecho(inferior_dir, direita, secao_rigida, ficticia=True)
            trecho(direita, superior_dir, secao_rigida, ficticia=True)
            atual = direita
        fim_viga = no(self.ap1 / 2.0 + self.l, 0.0)
        apoio_dir = no(self.ap1 / 2.0 + self.l + self.ap2 / 2.0, 0.0)
        apoio_dir.Tx = True
        apoio_dir.Ty = True
        trecho(atual, fim_viga, secao_cheia, carregado=True)
        trecho(fim_viga, apoio_dir, secao_cheia, carregado=True)

        self.portico = Estrutura(self.nos, self.barras)
        p = self.portico
        p.monta_k()
        p.monta_k01()
        p.monta_fnos()
        p.aplica_cc_fnos()
        p.calcula_deslocamentos()
        p.calcula_solicitacoes_internas_nodais()
        p.calcula_reacoes()
        for barra in self.barras:
            if not barra.ficticia:
                barra.calcula_vetor_ASL()
                barra.calcula_vetor_ASW_S()
        if imprimir_resultados:
            print('Deslocamentos:')
            print_vetor(p.u)
            print('Reações:')
            print_vetor(p.R)

    def visualizar(self, **kwargs):
        """Visualização de barras/diagramas (geometria de múltiplos furos requer adaptação)."""
        from .visualizacao import Visualizacao
        return Visualizacao(self).desenhar_geometria(**kwargs)
