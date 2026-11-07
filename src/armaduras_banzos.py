"""Pré-detalhamento longitudinal dos banzos dos furos (cm, cm²).

Uma proposta por face de cada banzo, usando a maior demanda ASL/ASL2
nas extremidades de seus elementos. Comprimentos geométricos apenas;
não incluem ancoragem, emendas ou prolongamentos de barras.
"""
from dataclasses import dataclass
from math import pi, isfinite
from .armaduras_longitudinais import _selecionar, BITOLAS_MM, propor_longitudinais


@dataclass(frozen=True)
class ArmaduraBanzo:
    posicao: str
    indice_furo: int
    banzo: str
    face: str
    as_necessaria_cm2: float
    as_fornecida_cm2: float
    quantidade: int
    bitola_mm: float
    inicio_cm: float
    fim_cm: float
    y_cm: float
    comprimento_geometrico_cm: float
    porta_estribo: bool = False
    as_total_necessaria_cm2: float = 0.0
    as_continua_cm2: float = 0.0

    @property
    def descricao(self):
        return (f'{self.posicao}: {self.quantidade} Ø{self.bitola_mm:g} '
                f'C={self.comprimento_geometrico_cm:.1f} cm '
                f'(As={self.as_fornecida_cm2:.2f} cm²)')


def _zonas(modelo):
    x0 = float(modelo.ap1) / 2
    h = float(modelo.h)
    if hasattr(modelo, 'furos'):
        furos = [(float(f.inicio), float(f.fim), float(f.a), float(f.c))
                 for f in modelo.furos]
    else:
        centro = float(modelo.lf)
        b = float(modelo.b)
        furos = [(centro - b / 2, centro + b / 2,
                  float(modelo.a), float(modelo.c))]
    for i, (inicio, fim, a, c) in enumerate(furos, 1):
        # Os banzos do modelo incluem as zonas rígidas de transição bw/2.
        esquerda = x0 + inicio - float(modelo.bw) / 2
        direita = x0 + fim + float(modelo.bw) / 2
        yield i, esquerda, direita, (-(h - c) / 2, c), ((a + c) / 2, h - a - c)


def propor_banzos(modelo, bitolas=BITOLAS_MM,
                   bitola_porta_estribo_mm=8.0, quantidade_porta_estribo=2,
                   primeira_posicao=3, cobrimento_esquematico_cm=3.5):
    """Escolhe barras complementares nos banzos, preservando N1/N2.

    Convenção idêntica ao detalhamento longitudinal original:
    momento positivo -> face inferior do próprio banzo; negativo -> superior.
    ASL2 fica na face oposta. A armadura contínua conta apenas nas faces
    externas dos banzos; faces internas sem demanda recebem porta-estribo.
    """
    # N1/N2 atravessam o furo: sua área fornecida contribui nas faces
    # EXTERNAS dos banzos (inferior do banzo inferior e superior do superior).
    continuas = {p.face: p.as_fornecida_cm2 for p in propor_longitudinais(modelo)}
    saida = []
    numero = primeira_posicao
    for indice, xa, xb, inferior, superior in _zonas(modelo):
        for nome, (yc, altura) in (('inferior', inferior), ('superior', superior)):
            demanda = {'inferior': 0.0, 'superior': 0.0}
            encontrados = 0
            for barra in modelo.barras:
                if getattr(barra, 'ficticia', False):
                    continue
                ni, nj = barra.noi, barra.noj
                x1, x2 = float(ni.x), float(nj.x)
                y1, y2 = float(ni.y), float(nj.y)
                if (abs(y1 - yc) > 1e-5 or abs(y2 - yc) > 1e-5
                        or abs(x2 - x1) < 1e-8
                        or min(x1, x2) < xa - 1e-5
                        or max(x1, x2) > xb + 1e-5):
                    continue
                momentos = (-float(barra.fl[2]), float(barra.fl[5]))
                for k, m in enumerate(momentos):
                    lado = 'inferior' if m >= 0 else 'superior'
                    oposto = 'superior' if lado == 'inferior' else 'inferior'
                    for face, as_req in ((lado, float(barra.ASL[k])),
                                         (oposto, float(barra.ASL2[k]))):
                        if not isfinite(as_req) or as_req < -1e-9:
                            raise ValueError(f'As inválida na barra {barra.num}.')
                        demanda[face] = max(demanda[face], max(0.0, as_req))
                encontrados += 1
            if not encontrados:
                raise ValueError(f'Nenhum elemento encontrado no banzo {nome} do furo {indice}.')
            for face in ('inferior', 'superior'):
                total_req = demanda[face]
                externa = ((nome == 'inferior' and face == 'inferior') or
                           (nome == 'superior' and face == 'superior'))
                existente = continuas.get(face, 0.0) if externa else 0.0
                area = max(0.0, total_req - existente)
                # Não criar barra adicional onde a armadura contínua já
                # atende à demanda. Em faces internas sem demanda, manter
                # o porta-estribo construtivo adotado anteriormente.
                if area <= 1e-10 and externa:
                    continue
                porta = area <= 1e-10
                if porta:
                    n = int(quantidade_porta_estribo)
                    diam = float(bitola_porta_estribo_mm)
                    if n < 2 or diam <= 0:
                        raise ValueError('Porta-estribo inválido.')
                    fornecida = n * pi * (diam / 10)**2 / 4
                else:
                    n, diam, fornecida = _selecionar(area, bitolas)
                deslocamento = min(cobrimento_esquematico_cm, altura * 0.25)
                y = yc + (altura / 2 - deslocamento if face == 'superior'
                          else -altura / 2 + deslocamento)
                saida.append(ArmaduraBanzo(
                    f'N{numero}', indice, nome, face, area, fornecida,
                    n, diam, xa, xb, y, xb - xa, porta,
                    total_req, existente))
                numero += 1
    return saida
