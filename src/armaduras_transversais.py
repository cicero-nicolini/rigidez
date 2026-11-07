"""Pré-detalhamento de estribos em trechos de seção cheia do rigidez.

Unidades: cm, cm²/m. Não detalha os banzos nem a região dos furos.
As propostas são preliminares: espaçamentos máximos normativos, ganchos,
ancoragem e acomodação devem ser verificados antes de uso executivo.
"""
from dataclasses import dataclass
from math import ceil, floor, isfinite, pi

BITOLAS_MM = (5.0, 6.3, 8.0, 10.0, 12.5)
ESPACAMENTOS_CM = (5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 30.0)


@dataclass(frozen=True)
class TrechoEstribos:
    posicao: str
    inicio_cm: float
    fim_cm: float
    asw_necessaria_cm2_m: float
    asw_fornecida_cm2_m: float
    bitola_mm: float
    ramos: int
    espacamento_cm: float
    quantidade: int

    @property
    def descricao(self):
        return (f'{self.posicao}: {self.quantidade} Ø{self.bitola_mm:g} '
                f'c/{self.espacamento_cm:g} cm '
                f'(Asw/s={self.asw_fornecida_cm2_m:.2f} cm²/m)')


def _intervalos_cheios(modelo):
    """Trechos de seção cheia, entre as faces dos apoios e os limites dos furos."""
    x0 = float(modelo.ap1) / 2
    x1 = x0 + float(modelo.l)
    intervalos = []
    atual = x0
    if hasattr(modelo, 'furos'):
        zonas = [(x0 + float(f.inicio) - float(modelo.bw) / 2,
                  x0 + float(f.fim) + float(modelo.bw) / 2)
                 for f in modelo.furos]
    elif all(hasattr(modelo, attr) for attr in ('lf', 'b', 'bw')):
        # Modelo legado (src.modelo.Modelo): lf medido da face do apoio.
        centro = x0 + float(modelo.lf)
        meia_zona = (float(modelo.b) + float(modelo.bw)) / 2
        zonas = [(centro - meia_zona, centro + meia_zona)]
    else:
        raise ValueError('Geometria das aberturas não identificada.')
    for a, b in sorted(zonas):
        if a > atual + 1e-7:
            intervalos.append((atual, a))
        atual = max(atual, b)
    if atual < x1 - 1e-7:
        intervalos.append((atual, x1))
    return intervalos


def _escolher(asw, bitolas, espacamentos, ramos, espacamento_maximo_cm):
    opcoes = []
    for d in bitolas:
        area = ramos * pi * (d / 10)**2 / 4
        for s in espacamentos:
            if s > espacamento_maximo_cm + 1e-8:
                continue
            fornecida = 100 * area / s
            if fornecida + 1e-9 >= asw:
                # Prefere menor excesso de aço e, em empate, maior espaçamento.
                opcoes.append((fornecida - asw, -s, d, s, fornecida))
    if not opcoes:
        raise ValueError(f'Não há estribo disponível para Asw/s={asw:.3f} cm²/m.')
    _, _, d, s, fornecida = min(opcoes)
    return d, s, fornecida


def propor_transversais(modelo, bitolas=BITOLAS_MM,
                        espacamentos=ESPACAMENTOS_CM, ramos=2,
                        espacamento_maximo_cm=20.0):
    """Uma proposta uniforme por trecho cheio, baseada no maior ASW_S local.

    Para manter a demanda conservadora nos limites, considera as extremidades
    dos elementos horizontais do eixo y=0 que pertencem ao intervalo.
    A primeira e última posições são preliminarmente colocadas nos limites.
    """
    if ramos < 2 or espacamento_maximo_cm <= 0:
        raise ValueError('Configuração de estribos inválida.')
    saida = []
    for indice, (a, b) in enumerate(_intervalos_cheios(modelo), 1):
        demanda = 0.0
        encontrados = 0
        for barra in modelo.barras:
            if getattr(barra, 'ficticia', False):
                continue
            ni, nj = barra.noi, barra.noj
            if abs(float(ni.y)) > 1e-6 or abs(float(nj.y)) > 1e-6:
                continue
            xa, xb = sorted((float(ni.x), float(nj.x)))
            if xb - xa < 1e-8 or xa < a - 1e-6 or xb > b + 1e-6:
                continue
            valores = list(barra.ASW_S)
            if len(valores) != 2 or any(not isfinite(float(v)) or float(v) < 0 for v in valores):
                raise ValueError(f'Barra {barra.num}: ASW_S inválido.')
            demanda = max(demanda, *(float(v) for v in valores))
            encontrados += 1
        if not encontrados:
            raise ValueError(f'Não foram encontrados elementos de seção cheia no trecho {indice}.')
        d, s, fornecida = _escolher(demanda, bitolas, espacamentos, ramos,
                                    espacamento_maximo_cm)
        # Quantidade para que a distância real não exceda o espaçamento escolhido.
        quantidade = max(2, ceil((b - a) / s - 1e-9) + 1)
        saida.append(TrechoEstribos(f'E{indice}', a, b, demanda, fornecida,
                                    d, ramos, s, quantidade))
    return saida
