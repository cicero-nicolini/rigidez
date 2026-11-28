"""Coordenação do detalhamento das armaduras de vigas."""

from .armaduras_longitudinais import propor_longitudinais
from .armaduras_banzos import propor_banzos
from .armaduras_transversais import propor_transversais
from .armaduras_estribos_banzos import propor_estribos_banzos
from .escalonamento_n1_n2 import propor_escalonamento, continuidade_nos_furos
from .ancoragem_apoios import ancorar_nos_apoios
from .ancoragem_banzos import ancorar_banzos
from .modulacao_armaduras import modular_armaduras

class Detalhamento:
    """Representa os resultados do detalhamento de uma viga."""

    def __init__(self, modelo) -> None:
        """
        Instancia um objeto da classe Detalhamento.

        Args:

            modelo: modelo de viga dimensionado
        """

        self.modelo = modelo
        self.armaduras_longitudinais = []
        self.armaduras_escalonadas = []
        self.armaduras_banzos = []
        self.armaduras_transversais = []
        self.armaduras_estribos_banzos = []
        self.diagnosticos_ancoragem_apoios = []
        self.diagnosticos_ancoragem_banzos = []

    def calcula_longitudinais(self, niveis_n1:int = 0, niveis_n2:int = 0,
                              decalagem_cm:float = None, decalagem_automatica:bool = True,
                              ancoragem_automatica:bool = False, ancoragem_apoios:bool = True,
                              boa_aderencia:bool = True) -> None:
        """
        Calcula as armaduras longitudinais, escalonamento e ancoragens.

        Args:

            niveis_n1: níveis adicionais inferiores
            niveis_n2: níveis adicionais superiores
            decalagem_cm: distância manual em cm
            decalagem_automatica: habilita decalagem automática
            ancoragem_automatica: habilita cálculo de ancoragem
            ancoragem_apoios: habilita ancoragem nos apoios
            boa_aderencia: indica condição de boa aderência
        """

        m = self.modelo
        limite_x_esquerdo = float(m.ap1)/2.0
        limite_x_direito = limite_x_esquerdo + float(m.l)
        self.armaduras_longitudinais = propor_longitudinais(m)
        self.diagnosticos_ancoragem_banzos = []

        faixas = propor_escalonamento(
            m, niveis_n1=niveis_n1, niveis_n2=niveis_n2,
            decalagem_cm=decalagem_cm,
            decalagem_automatica=decalagem_automatica,
            ancoragem_automatica=ancoragem_automatica,
            boa_aderencia=boa_aderencia)
        faixas, areas_continuas = continuidade_nos_furos(m, faixas)
        if ancoragem_automatica and ancoragem_apoios:

            faixas, self.diagnosticos_ancoragem_apoios = ancorar_nos_apoios(
                m, faixas, boa_aderencia=boa_aderencia)
        else:
            self.diagnosticos_ancoragem_apoios = []
        # Permite prolongamento até as faces externas dos apoios.
        # Se a modulação não couber, a rotina sinaliza insuficiência.
        # As barras retas podem entrar nos apoios, mas devem respeitar
        # o cobrimento na face externa de cada pilar.
        cobrimento_apoio = float(getattr(m, 'cob', 0.0))
        # O vão l é livre entre faces internas dos pilares.
        # A ancoragem parte dessas faces, não dos nós centrais.
        limite_mod_esq = (limite_x_esquerdo - float(m.ap1)
                          + cobrimento_apoio
                          if ancoragem_automatica and ancoragem_apoios
                          else limite_x_esquerdo)
        limite_mod_dir = (limite_x_direito + float(m.ap2)
                          - cobrimento_apoio
                          if ancoragem_automatica and ancoragem_apoios
                          else limite_x_direito)
        faixas = modular_armaduras(
            faixas, limite_mod_esq, limite_mod_dir)
        self.armaduras_escalonadas = faixas
        # Redimensiona apenas os complementos dos banzos externos,
        # descontando a armadura N1/N2 que atravessa cada furo.
        novos_banzos = propor_banzos(
            m, areas_continuas_por_furo=areas_continuas)
        if ancoragem_automatica:

            novos_banzos, self.diagnosticos_ancoragem_banzos = ancorar_banzos(
                m, novos_banzos, boa_aderencia=boa_aderencia)
        if ancoragem_automatica and ancoragem_apoios:
            novos_banzos, diagnosticos_banzos_apoios = ancorar_nos_apoios(
                m, novos_banzos, boa_aderencia=boa_aderencia)
            self.diagnosticos_ancoragem_apoios.extend(diagnosticos_banzos_apoios)
        novos_banzos = modular_armaduras(
            novos_banzos, limite_mod_esq, limite_mod_dir)
        self.armaduras_banzos = novos_banzos


    def calcula_transversais(self) -> None:
        """Calcula as armaduras transversais da seção cheia e dos banzos."""

        m = self.modelo
        self.armaduras_transversais = propor_transversais(m)
        self.armaduras_estribos_banzos = propor_estribos_banzos(
            m, primeira_posicao=len(self.armaduras_transversais) + 1)
