"""Validação de ações uniformemente distribuídas e combinações simples."""
import math


def normalizar_cargas(cargas=None, w=None):
    if cargas is None:
        if w is None:
            raise ValueError("Informe cargas ou w (compatibilidade).")
        cargas = [{"nome": "Carga distribuída", "w": w, "gamaf": 1.4}]
    if not isinstance(cargas, (list, tuple)) or not cargas:
        raise ValueError("cargas deve ser uma lista não vazia.")
    resultado = []
    for i, item in enumerate(cargas, 1):
        if not isinstance(item, dict):
            raise TypeError(f"Carga {i}: esperado dicionário.")
        nome = str(item.get("nome", f"Carga {i}"))
        try:
            intensidade = float(item["w"])
            fator = float(item["gamaf"])
        except (KeyError, ValueError, TypeError) as exc:
            raise ValueError(f"Carga {i}: informe w e gamaf numéricos.") from exc
        if not math.isfinite(intensidade) or not math.isfinite(fator) or fator < 0:
            raise ValueError(f"Carga {i}: w e gamaf devem ser finitos; gamaf >= 0.")
        resultado.append({"nome": nome, "w": intensidade, "gamaf": fator})
    return resultado


def analisar_combinacoes(modelo, imprimir=False):
    """Calcula dois estados lineares, deixando o estado de cálculo ativo.

    O estado característico é armazenado em resultados_caracteristicos;
    o de cálculo em resultados_calculo. Dimensionamento usa gamaf=1.
    """
    import numpy as np
    from .carregamento_distribuido import Carregamento_distribuido

    p = modelo.portico
    modelo.resultados_caracteristicos = {}
    modelo.resultados_calculo = {}
    for chave, intensidade in (("caracteristicos", modelo.wk), ("calculo", modelo.wd)):
        for barra in modelo.barras:
            barra.fepl[:] = 0
            barra.fep[:] = 0
        modelo.carregamentos = []
        for barra in modelo.barras:
            if getattr(barra, "carregada", False):
                carga = Carregamento_distribuido(0, barra.L, intensidade, intensidade, barra)
                modelo.carregamentos.append(carga)
                carga.calcula_fepl()
                carga.calcula_fep()
        p.monta_fnos()
        p.aplica_cc_fnos()
        p.calcula_deslocamentos()
        p.calcula_solicitacoes_internas_nodais()
        p.calcula_reacoes()
        resultado = {
            "u": np.array(p.u, copy=True),
            "R": np.array(p.R, copy=True),
            "fl": {barra.num: barra.fl.copy() for barra in modelo.barras},
            "f": {barra.num: barra.f.copy() for barra in modelo.barras},
            "fep": {barra.num: barra.fep.copy() for barra in modelo.barras},
        }
        setattr(modelo, "resultados_" + chave, resultado)
    modelo.tipo_esforcos = "calculo"
    for barra in modelo.barras:
        if not barra.ficticia:
            barra.gamaf = 1.0  # fl já contém os esforços majorados
            barra.calcula_vetor_ASL()
            barra.calcula_vetor_ASW_S()


def selecionar_esforcos(modelo, tipo="calculo"):
    """Altera o estado exibido (fl), sem recalcular armaduras.

    Para voltar ao dimensionamento use tipo='calculo'.
    """
    if tipo not in ("calculo", "caracteristicos"):
        raise ValueError("tipo deve ser 'calculo' ou 'caracteristicos'.")
    estado = getattr(modelo, "resultados_" + tipo)
    modelo.portico.u = estado["u"].copy()
    modelo.portico.R = estado["R"].copy()
    for barra in modelo.barras:
        barra.fl = estado["fl"][barra.num].copy()
        barra.f = estado["f"][barra.num].copy()
    modelo.tipo_esforcos = tipo
