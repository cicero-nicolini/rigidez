"""Pré-detalhamento de ancoragem reta das armaduras longitudinais dos banzos.

Calcula l_b,nec e prolonga geometricamente as barras além das faces dos furos.
Não verifica disponibilidade física de ancoragem, confinamento, emendas,
interferências nem regras específicas de ancoragem em apoios.
"""
from dataclasses import dataclass, replace
from math import isfinite

from .escalonamento_n1_n2 import comprimento_ancoragem_necessario
from .armaduras_banzos import _zonas


@dataclass(frozen=True)
class DiagnosticoAncoragemBanzo:
    posicao: str
    lb_necessario_cm: float
    esquerda_disponivel_cm: float
    direita_disponivel_cm: float
    esquerda_atendida: bool
    direita_atendida: bool


def ancorar_banzos(modelo, propostas, boa_aderencia=True):
    """Prolonga cada armadura de banzo por l_b,nec em ambas as extremidades.

    Args:
        modelo: Modelo estrutural com barras, materiais e geometria em cm.
        propostas: Propostas de ``propor_banzos``.
        boa_aderencia: Condição de aderência para o cálculo de l_b,nec.

    Returns:
        Tupla (propostas_atualizadas, diagnosticos). As áreas não são alteradas.

    Notes:
        Para porta-estribos sem demanda de cálculo, adota-se As,calc/As,ef=1
        como critério conservador preliminar. As extensões são limitadas ao
        domínio longitudinal do modelo; uma extensão truncada NÃO constitui
        ancoragem verificada.
    """
    x_min = float(modelo.ap1) / 2.0
    x_max = x_min + float(modelo.l)
    zonas = {indice: {'inferior': inferior[0], 'superior': superior[0]}
             for indice, _, _, inferior, superior in _zonas(modelo)}
    novas, diagnosticos = [], []
    for p in propostas:
        y_centro = zonas[p.indice_furo][p.banzo]
        candidatos = []
        for barra in modelo.barras:
            if getattr(barra, 'ficticia', False):
                continue
            no1, no2 = barra.noi, barra.noj
            if (abs(float(no1.y) - float(no2.y)) > 1e-5 or
                    abs(float(no1.y) - y_centro) > 1e-5):
                continue
            xa, xb = sorted((float(no1.x), float(no2.x)))
            if xa >= float(p.inicio_cm) - 1e-5 and xb <= float(p.fim_cm) + 1e-5:
                candidatos.append(barra)
        if not candidatos:
            raise ValueError(f'Nenhuma barra estrutural associada a {p.posicao}.')
        # Maior comprimento exigido entre os elementos do banzo.
        as_calc = (float(p.as_necessaria_cm2) if not p.porta_estribo
                   else float(p.as_fornecida_cm2))
        lb = max(comprimento_ancoragem_necessario(
            barra, p.bitola_mm, as_calc, p.as_fornecida_cm2,
            boa_aderencia=boa_aderencia) for barra in candidatos)
        inicio = max(x_min, float(p.inicio_cm) - lb)
        fim = min(x_max, float(p.fim_cm) + lb)
        disp_esq = float(p.inicio_cm) - inicio
        disp_dir = fim - float(p.fim_cm)
        novas.append(replace(p, inicio_cm=inicio, fim_cm=fim,
                             comprimento_geometrico_cm=fim - inicio))
        diagnosticos.append(DiagnosticoAncoragemBanzo(
            p.posicao, lb, disp_esq, disp_dir,
            disp_esq >= lb - 1e-6, disp_dir >= lb - 1e-6))
    return novas, diagnosticos
