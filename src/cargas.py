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


def analisar_rigidez_banzos(modelo, max_iter:int = 30, tolerancia:float = 1e-8) -> int:
    """
    Atualiza iterativamente a rigidez dos banzos pelo esforço normal médio.

    Args:

        modelo: modelo com carregamentos da combinação atual
        max_iter: número máximo de análises
        tolerancia: tolerância absoluta do esforço normal (kN)

    Returns:

        Número de análises realizadas até estabilizar os modificadores.

    Raises:

        RuntimeError: ciclo de estados ou ausência de convergência.
    """
    p = modelo.portico
    grupos = getattr(modelo, 'banzos_furos', ())
    for barra in modelo.barras:
        barra.beta_flexao = 1.0

    visitados = set()
    for iteracao in range(1, max_iter + 1):
        atual = tuple(b.beta_flexao for b in modelo.barras)
        if atual in visitados:
            raise RuntimeError('Oscilação na identificação dos banzos tracionados.')
        visitados.add(atual)

        p.monta_k()
        p.monta_k01()
        p.monta_fnos()
        p.aplica_cc_fnos()
        p.calcula_deslocamentos()
        p.calcula_solicitacoes_internas_nodais()

        # Normal positivo em tração: N inicial = -fl[0], N final = fl[3].
        novos = {b.num: 1.0 for b in modelo.barras}
        for furo, inferior, superior in grupos:
            for banzo in (inferior, superior):
                normal = sum(b.L*(b.fl[3]-b.fl[0])/2 for b in banzo)/sum(b.L for b in banzo)
                if normal > tolerancia:
                    for b in banzo:
                        novos[b.num] = furo.beta_t

        if atual == tuple(novos[b.num] for b in modelo.barras):
            p.calcula_reacoes()
            return iteracao

        for barra in modelo.barras:
            barra.beta_flexao = novos[barra.num]

    raise RuntimeError('Rigidez dos banzos não convergiu no limite de iterações.')


def analisar_combinacoes(modelo, imprimir=False):
    """Calcula dois estados de carregamento, deixando o estado de cálculo ativo.

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
        iteracoes = analisar_rigidez_banzos(modelo)
        resultado = {
            "iteracoes": iteracoes,
            "beta_flexao": {b.num: b.beta_flexao for b in modelo.barras},
            "carregamentos": tuple(modelo.carregamentos),
            "u": np.array(p.u, copy=True),
            "R": np.array(p.R, copy=True),
            "fnos": p.fnos.copy(),
            "f_global": p.f.copy(),
            "u_barras": {barra.num: barra.u.copy() for barra in modelo.barras},
            "fepl": {barra.num: barra.fepl.copy() for barra in modelo.barras},
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
    """
    Restaura os resultados da combinação selecionada sem recalcular armaduras.

    Args:

        modelo: objeto modelo previamente analisado
        tipo: combinação de cálculo ou característica

    Raises:

        ValueError: tipo de combinação inválido.
    """
    if tipo not in ("calculo", "caracteristicos"):
        raise ValueError("tipo deve ser 'calculo' ou 'caracteristicos'.")
    estado = getattr(modelo, "resultados_" + tipo)
    modelo.carregamentos = list(estado["carregamentos"])
    modelo.portico.u = estado["u"].copy()
    modelo.portico.R = estado["R"].copy()
    modelo.portico.fnos = estado["fnos"].copy()
    modelo.portico.f = estado["f_global"].copy()

    # Restaura os vetores de cada barra no estado selecionado
    for barra in modelo.barras:
        barra.fl = estado["fl"][barra.num].copy()
        barra.f = estado["f"][barra.num].copy()
        barra.u = estado["u_barras"][barra.num].copy()
        barra.fep = estado["fep"][barra.num].copy()
        barra.fepl = estado["fepl"][barra.num].copy()
    # Restaura também as matrizes associadas à rigidez da combinação
    for barra in modelo.barras:
        barra.beta_flexao = estado["beta_flexao"][barra.num]
    modelo.portico.monta_k()
    modelo.portico.monta_k01()
    modelo.tipo_esforcos = tipo
