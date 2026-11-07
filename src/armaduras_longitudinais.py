"""Pré-detalhamento longitudinal do modelo rigidez (unidades: cm e cm²).

As barras propostas cobrem todo o vão entre faces internas dos apoios.
Comprimento geométrico NÃO inclui ancoragens, dobras ou emendas.
Apenas trechos horizontais de seção cheia são usados para obter As,req.
"""
from dataclasses import dataclass
from math import ceil, pi, isfinite

BITOLAS_MM = (8.0, 10.0, 12.5, 16.0, 20.0, 25.0, 32.0)


@dataclass(frozen=True)
class ArmaduraLongitudinal:
    posicao: str
    face: str
    as_necessaria_cm2: float
    as_fornecida_cm2: float
    quantidade: int
    bitola_mm: float
    inicio_cm: float
    fim_cm: float
    comprimento_geometrico_cm: float

    @property
    def descricao(self):
        return (f'{self.posicao}: {self.quantidade} Ø{self.bitola_mm:g} '
                f'C={self.comprimento_geometrico_cm:.1f} cm '
                f'(As={self.as_fornecida_cm2:.2f} cm²)')


def _selecionar(as_req, bitolas=BITOLAS_MM, minimo=2, maximo=12):
    """Menor excesso de aço; empate favorece menos barras e menor bitola."""
    opcoes = []
    for diam in bitolas:
        area = pi * (diam / 10.0)**2 / 4.0
        n = max(minimo, ceil((as_req - 1e-10) / area))
        if n <= maximo:
            opcoes.append((n * area - as_req, n, diam, n * area))
    if not opcoes:
        raise ValueError('Área necessária excede as opções automáticas; revisar bitolas/camadas.')
    _, n, diam, fornecida = min(opcoes)
    return n, diam, fornecida


def propor_longitudinais(modelo, bitolas=BITOLAS_MM,
                         bitola_porta_estribo_mm=8.0, quantidade_porta_estribo=2):
    """Propõe uma armadura contínua por face para a demanda máxima em seção cheia.

    A convenção de sinais segue src.visualizacao: M=-fl[2] na origem e
    M=fl[5] no fim; M>0 -> face inferior; M<0 -> face superior.
    ASL2, quando presente, fica na face oposta à ASL.
    """
    x0 = float(modelo.ap1) / 2.0
    x1 = x0 + float(modelo.l)
    if not (isfinite(x0) and isfinite(x1) and x1 > x0):
        raise ValueError('Vão inválido.')
    demanda = {'inferior': 0.0, 'superior': 0.0}
    selecionadas = 0
    for barra in modelo.barras:
        if getattr(barra, 'ficticia', False):
            continue
        ni, nj = barra.noi, barra.noj
        xa, xb = float(ni.x), float(nj.x)
        # Apenas barras no eixo da seção cheia, dentro do vão.
        if abs(float(ni.y)) > 1e-6 or abs(float(nj.y)) > 1e-6:
            continue
        if min(xa, xb) < x0 - 1e-6 or max(xa, xb) > x1 + 1e-6:
            continue
        if abs(xb - xa) < 1e-8:
            continue
        fl = barra.fl
        if len(fl) < 6:
            raise ValueError(f'Barra {barra.num}: fl inválido.')
        momentos = (-float(fl[2]), float(fl[5]))
        as1, as2 = barra.ASL, barra.ASL2
        for k, m in enumerate(momentos):
            if not isfinite(m):
                raise ValueError(f'Barra {barra.num}: momento não finito.')
            lado = 'inferior' if m >= 0 else 'superior'
            oposto = 'superior' if lado == 'inferior' else 'inferior'
            for face, valor in ((lado, float(as1[k])), (oposto, float(as2[k]))):
                if not isfinite(valor) or valor < -1e-9:
                    raise ValueError(f'Barra {barra.num}: área de aço inválida.')
                demanda[face] = max(demanda[face], max(0.0, valor))
        selecionadas += 1
    if not selecionadas:
        raise ValueError('Não foram encontradas barras de seção cheia dentro do vão.')
    saida = []
    for posicao, face in (('N1', 'inferior'), ('N2', 'superior')):
        area = demanda[face]
        if area <= 1e-10:
            if face != 'superior':
                continue
            # Armadura construtiva porta-estribos, quando não há demanda
            # calculada na face superior. Bitola futuramente vinculada
            # à seleção dos estribos.
            n = int(quantidade_porta_estribo)
            diam = float(bitola_porta_estribo_mm)
            if n < 2 or not isfinite(diam) or diam <= 0:
                raise ValueError('Configuração inválida do porta-estribos.')
            fornecida = n * pi * (diam / 10.0)**2 / 4.0
        else:
            n, diam, fornecida = _selecionar(area, bitolas)
        saida.append(ArmaduraLongitudinal(
            posicao, face, area, fornecida, n, diam,
            x0, x1, x1 - x0))
    return saida
