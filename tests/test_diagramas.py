"""Validação do pós-processamento de deslocamentos e esforços."""

import numpy as np
import pytest
from conftest import analisar
from src.carregamento_pontual import Carregamento_pontual


def test_deformada_balanco(balanco):
    """Compara a interpolação cúbica com a solução de Euler-Bernoulli."""
    p = balanco
    b = p.barras[0]
    F = 20.
    p.nos[1].Fy = -F
    analisar(p)
    _, xy, u = p.calcula_deformada(21)[0]
    x = xy[:, 0]
    esperado = -F*x**2*(3*b.L-x)/(6*b.secao.E*b.secao.I)
    np.testing.assert_allclose(u[:, 1], esperado, atol=1e-10)
    np.testing.assert_allclose(u[:, 0], 0., atol=1e-10)


def test_diagrama_uniforme(biapoiada):
    """Compara N, V e M com o equilíbrio analítico da viga."""
    p, cargas = biapoiada
    q, L = cargas[0].w1, 600.
    for b, x, N, V, M in p.calcula_solicitacoes_internas(cargas):
        X = x+b.noi.x
        np.testing.assert_allclose(N, 0., atol=1e-8)
        np.testing.assert_allclose(V, q*(X-L/2), atol=1e-8)
        np.testing.assert_allclose(M, q*X*(L-X)/2, atol=1e-8)


def test_salto_cortante(balanco):
    """Verifica a descontinuidade do cortante sob força concentrada."""
    p = balanco
    b = p.barras[0]
    carga = Carregamento_pontual(120., 0., 20., b)
    carga.calcula_fepl()
    carga.calcula_fep()
    analisar(p)
    _, x, N, V, M = p.calcula_solicitacoes_internas([carga])[0]
    indices = np.flatnonzero(np.isclose(x, 120.))
    assert len(indices) == 2
    assert V[indices[1]]-V[indices[0]] == pytest.approx(20.)
    assert M[indices[1]] == pytest.approx(M[indices[0]], abs=1e-8)


def test_deformada_pontos_invalidos(balanco):
    """Rejeita discretização com menos de dois pontos."""
    with pytest.raises(ValueError):
        balanco.calcula_deformada(1)
