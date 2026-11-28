"""Validação da análise matricial por soluções analíticas e regressão."""

import numpy as np
import pytest
from conftest import analisar


def test_tracao(balanco):
    """Verifica deslocamento axial e reação de uma barra tracionada."""
    p = balanco
    b = p.barras[0]
    p.nos[1].Fx = 50.
    analisar(p)
    assert p.u[3] == pytest.approx(50*b.L/(b.secao.E*b.secao.A))
    assert p.R[0] == pytest.approx(-50.)
    np.testing.assert_allclose(b.fl, [-50, 0, 0, 50, 0, 0], atol=1e-8)


def test_balanco(balanco):
    """Verifica flecha, rotação e reações de uma viga em balanço."""
    p = balanco
    b = p.barras[0]
    F, L, EI = 20., b.L, b.secao.E*b.secao.I
    p.nos[1].Fy = -F
    analisar(p)
    assert p.u[4] == pytest.approx(-F*L**3/(3*EI))
    assert p.u[5] == pytest.approx(-F*L**2/(2*EI))
    assert p.R[1] == pytest.approx(F)
    assert p.R[2] == pytest.approx(F*L)


def test_biapoiada(biapoiada):
    """Verifica flecha central e reações sob carga uniforme."""
    p, cargas = biapoiada
    L, q = 600., cargas[0].w1
    s = p.barras[0].secao
    assert p.u[4] == pytest.approx(-5*q*L**4/(384*s.E*s.I))
    np.testing.assert_allclose(p.R[[1, 7]], [q*L/2]*2, atol=1e-8)
    assert p.R[1]+p.R[7] == pytest.approx(q*L)


def test_matriz_original(balanco):
    """Verifica simetria, independência e preservação da matriz global."""
    p = balanco
    original = p.k.copy()
    p.monta_k01()
    assert not np.shares_memory(p.k, p.k01)
    np.testing.assert_allclose(p.k, original)
    np.testing.assert_allclose(p.k, p.k.T)
    np.testing.assert_allclose(p.k01[:3], np.eye(6)[:3])
    np.testing.assert_allclose(p.k01[:, :3], np.eye(6)[:, :3])


def test_remontagem(balanco):
    """Verifica que a remontagem não acumula rigidez."""
    original = balanco.k.copy()
    for _ in range(3):
        balanco.monta_k()
        np.testing.assert_allclose(balanco.k, original)


def test_reacoes_repetidas(biapoiada):
    """Verifica que as reações não acumulam em chamadas sucessivas."""
    p, _ = biapoiada
    original = p.R.copy()
    for _ in range(3):
        p.calcula_reacoes()
        np.testing.assert_allclose(p.R, original, atol=1e-8)


def test_residuo(biapoiada):
    """Verifica o resíduo do sistema com condições de contorno."""
    p, _ = biapoiada
    np.testing.assert_allclose(p.k01 @ p.u, p.fnos, atol=1e-8)
