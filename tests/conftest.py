"""Configuração e modelos auxiliares para os testes numéricos."""

import pytest
from src.no import No
from src.barra import Barra
from src.secao import Secao
from src.estrutura import Estrutura
from src.carregamento_distribuido import Carregamento_distribuido
from src.furo import Furo
from src.modelo_furos import ModeloFuros


def analisar(portico: Estrutura) -> None:
    """Executa a análise de uma estrutura com rigidez e cargas já definidas."""

    portico.monta_fnos()
    portico.aplica_cc_fnos()
    portico.calcula_deslocamentos()
    portico.calcula_solicitacoes_internas_nodais()
    portico.calcula_reacoes()


@pytest.fixture
def balanco() -> Estrutura:
    """Cria uma barra horizontal engastada à esquerda (cm, kN)."""

    a, b = No(1, 0., 0.), No(2, 300., 0.)
    a.Tx = a.Ty = a.Rz = True
    secao = Secao('retangular', 20., 40., 3000.)
    p = Estrutura([a, b], [Barra(1, a, b, secao)])
    p.monta_k()
    p.monta_k01()
    return p


@pytest.fixture
def biapoiada() -> tuple:
    """Cria uma viga biapoiada com carga uniforme e dois elementos."""

    nos = [No(i+1, x, 0.) for i, x in enumerate((0., 300., 600.))]
    nos[0].Tx = nos[0].Ty = nos[-1].Ty = True
    secao = Secao('retangular', 20., 40., 3000.)
    barras = [Barra(i+1, nos[i], nos[i+1], secao) for i in range(2)]
    p = Estrutura(nos, barras)
    p.monta_k()
    p.monta_k01()
    cargas = [Carregamento_distribuido(0., b.L, .2, .2, b) for b in barras]
    for carga in cargas:
        carga.calcula_fepl()
        carga.calcula_fep()
    analisar(p)
    return p, cargas


@pytest.fixture(scope='module')
def modelo() -> ModeloFuros:
    """Cria uma viga com abertura e duas combinações de carregamento."""

    return ModeloFuros(bw=30., h=80., l=600.,
                       furos=[Furo(lf=250., b=100., a=30., c=25.)],
                       ap1=40., ap2=40., fck=35., fyk=500., fywk=500.,
                       cargas=[{'w': .15, 'gamaf': 1.4},
                               {'w': .10, 'gamaf': 1.5}],
                       espacamento_maximo=30.)
