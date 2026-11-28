"""Testes de regressão dos estados característico e de cálculo."""

import numpy as np
import pytest
from conftest import analisar
from src.cargas import normalizar_cargas


def test_combinacoes_lineares(modelo):
    """Verifica proporcionalidade dos resultados para cargas uniformes."""
    k, d = modelo.resultados_caracteristicos, modelo.resultados_calculo
    fator = modelo.wd/modelo.wk
    for chave in ('u', 'R', 'fnos'):
        np.testing.assert_allclose(d[chave], fator*k[chave], atol=1e-7)
    for b in modelo.barras:
        np.testing.assert_allclose(d['fl'][b.num], fator*k['fl'][b.num], atol=1e-7)


def test_restauracao(modelo):
    """Verifica a restauração completa e o recálculo de cada combinação."""
    p = modelo.portico
    armaduras = [(b.ASL.copy(), b.ASL2.copy(), b.ASW_S.copy()) for b in modelo.barras]
    try:
        for tipo in ('caracteristicos', 'calculo', 'caracteristicos'):
            modelo.selecionar_esforcos(tipo)
            estado = getattr(modelo, 'resultados_'+tipo)
            assert modelo.tipo_esforcos == tipo
            assert tuple(modelo.carregamentos) == estado['carregamentos']
            for atributo, chave in [('u','u'), ('R','R'), ('fnos','fnos'), ('f','f_global')]:
                atual = getattr(p, atributo)
                np.testing.assert_allclose(atual, estado[chave])
                assert not np.shares_memory(atual, estado[chave])
            for b, areas in zip(modelo.barras, armaduras):
                for atributo, chave in [('u','u_barras'), ('fl','fl'), ('f','f'), ('fep','fep'), ('fepl','fepl')]:
                    atual = getattr(b, atributo)
                    np.testing.assert_allclose(atual, estado[chave][b.num])
                    assert not np.shares_memory(atual, estado[chave][b.num])
                for atual, esperado in zip((b.ASL,b.ASL2,b.ASW_S), areas):
                    np.testing.assert_allclose(atual, esperado)
            analisar(p)
            np.testing.assert_allclose(p.R, estado['R'], atol=1e-7)
            np.testing.assert_allclose(p.u, estado['u'], atol=1e-9)
            for b in modelo.barras:
                np.testing.assert_allclose(b.fl, estado['fl'][b.num], atol=1e-7)
    finally:
        modelo.selecionar_esforcos('calculo')


def test_combinacao_invalida(modelo):
    """Rejeita um identificador de combinação desconhecido."""
    with pytest.raises(ValueError):
        modelo.selecionar_esforcos('invalido')


@pytest.mark.parametrize('cargas', [[], [{'w': 1, 'gamaf': -1}], [{'w': float('nan'), 'gamaf': 1}], [{'w': 1}]])
def test_cargas_invalidas(cargas):
    """Rejeita listas vazias, fatores negativos e dados não finitos ou ausentes."""
    with pytest.raises(ValueError):
        normalizar_cargas(cargas)
