"""Proposta preliminar de estribos dos banzos de furos (cm, cm²/m)."""
from dataclasses import dataclass
from math import ceil, isfinite
from .armaduras_banzos import _zonas
from .armaduras_transversais import _escolher, BITOLAS_MM, ESPACAMENTOS_CM

@dataclass(frozen=True)
class EstribosBanzo:
    posicao: str
    indice_furo: int
    banzo: str
    inicio_cm: float
    fim_cm: float
    y_inferior_cm: float
    y_superior_cm: float
    asw_necessaria_cm2_m: float
    asw_fornecida_cm2_m: float
    bitola_mm: float
    ramos: int
    espacamento_cm: float
    quantidade: int

    @property
    def descricao(self):
        return (f'{self.posicao} Ø{self.bitola_mm:g} c/{self.espacamento_cm:g} '
                f'(Asws={self.asw_fornecida_cm2_m:.2f} cm²/m)')


def propor_estribos_banzos(modelo, primeira_posicao=3,
                           bitolas=BITOLAS_MM, espacamentos=ESPACAMENTOS_CM,
                           ramos=2, espacamento_maximo_cm=20.0,
                           cobrimento_esquematico_cm=3.5):
    """Escolhe uma posição E# por banzo, com base no máximo ASW_S.

    Não verifica os critérios normativos de detalhamento executivo.
    """
    resultado = []
    numero = primeira_posicao
    for indice, xa, xb, inferior, superior in _zonas(modelo):
        for nome, (yc, altura) in (('inferior', inferior), ('superior', superior)):
            demanda = 0.0
            encontrados = 0
            for barra in modelo.barras:
                if getattr(barra, 'ficticia', False):
                    continue
                ni, nj = barra.noi, barra.noj
                if abs(float(ni.y) - yc) > 1e-5 or abs(float(nj.y) - yc) > 1e-5:
                    continue
                x1, x2 = sorted((float(ni.x), float(nj.x)))
                if x2-x1 < 1e-8 or x1 < xa-1e-5 or x2 > xb+1e-5:
                    continue
                valores = [float(v) for v in barra.ASW_S]
                if len(valores) != 2 or any(not isfinite(v) or v < 0 for v in valores):
                    raise ValueError(f'ASW_S inválido na barra {barra.num}.')
                demanda = max(demanda, *valores)
                encontrados += 1
            if not encontrados:
                raise ValueError(f'Sem elementos horizontais no banzo {nome} do furo {indice}.')
            d, s, fornecida = _escolher(demanda, bitolas, espacamentos,
                                        ramos, espacamento_maximo_cm)
            quantidade = max(2, ceil((xb-xa)/s - 1e-9) + 1)
            cobrimento = min(cobrimento_esquematico_cm, altura * .2)
            resultado.append(EstribosBanzo(
                f'E{numero}', indice, nome, xa, xb,
                yc-altura/2+cobrimento, yc+altura/2-cobrimento,
                demanda, fornecida, d, ramos, s, quantidade))
            numero += 1
    return resultado
