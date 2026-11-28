"""Testes dos modificadores de rigidez à flexão dos banzos."""
import numpy as np
import pytest
from src.furo import Furo
from src.modelo_furos import ModeloFuros


def criar_modelo(furos, **kwargs):
    dados = dict(bw=30., h=80., l=1100., furos=furos, ap1=40., ap2=40.,
                 w=.45, fck=35., fyk=500., fywk=500., espacamento_maximo=40.,
                 vinculacao_esquerda=1, vinculacao_direita=3)
    dados.update(kwargs)
    return ModeloFuros(**dados)


@pytest.mark.parametrize('beta', [0., -1., 1.1, float('nan')])
def test_beta_invalido(beta):
    with pytest.raises(ValueError):
        Furo(200, 100, 30, 25, beta)


def test_rigidez_elementar():
    m = criar_modelo([Furo(200, 100, 30, 25)])
    b = m.banzos_furos[0][1][0]
    original = b.kl.copy()
    b.beta_flexao = .4
    b.calcula_klocal()
    axial, flexao = [0, 3], [1, 2, 4, 5]
    np.testing.assert_allclose(b.kl[np.ix_(axial, axial)], original[np.ix_(axial, axial)])
    np.testing.assert_allclose(b.kl[np.ix_(flexao, flexao)], .4*original[np.ix_(flexao, flexao)])


def test_multiplos_furos():
    m = criar_modelo([Furo(200, 100, 30, 25, .5), Furo(500, 120, 25, 30, .3),
                      Furo(1000, 80, 35, 20, .7)])
    lados = set()
    for f, inf, sup in m.banzos_furos:
        for nome, banzo in [('inf', inf), ('sup', sup)]:
            normal = sum(b.L*(b.fl[3]-b.fl[0])/2 for b in banzo)/sum(b.L for b in banzo)
            esperado = f.beta_t if normal > 1e-8 else 1.
            assert all(b.beta_flexao == esperado for b in banzo)
            if esperado < 1:
                lados.add(nome)
    assert lados == {'inf', 'sup'}
    assert all(b.beta_flexao == 1 for b in m.barras if b.ficticia)
    assert m.resultados_calculo['iteracoes'] >= 2


def test_restauracao_rigidez():
    m = criar_modelo([Furo(200, 100, 30, 25, .5)],
                     cargas=[{'w': .1, 'gamaf': 1}, {'w': -.08, 'gamaf': 2}])
    assert m.resultados_calculo['beta_flexao'] != m.resultados_caracteristicos['beta_flexao']
    for tipo in ('caracteristicos', 'calculo'):
        m.selecionar_esforcos(tipo)
        p = m.portico
        estado = getattr(m, 'resultados_'+tipo)
        assert all(b.beta_flexao == estado['beta_flexao'][b.num] for b in m.barras)
        np.testing.assert_allclose(p.k01 @ p.u, p.fnos, atol=1e-7)
