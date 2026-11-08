"""Pré-detalhamento de barras retas nas regiões dos apoios.

Verifica ancoragem reta e propõe gancho com alpha=0,7 quando necessário,
sem presumir que
uma ancoragem insuficiente esteja atendida. NÃO dimensiona ganchos,
dobras, nós viga-pilar, confinamento nem ancoragens normativas específicas.
Unidades geométricas: cm.
"""
from dataclasses import dataclass, replace
from math import isfinite, pi

from .escalonamento_n1_n2 import comprimento_ancoragem_necessario


@dataclass(frozen=True)
class DiagnosticoAncoragemApoio:
    posicao: str
    lado: str
    face: str
    lb_necessario_cm: float
    comprimento_disponivel_cm: float
    comprimento_adotado_cm: float
    atendida_geometricamente: bool
    observacao: str
    tipo_ancoragem: str = "reta"
    lb_gancho_cm: float | None = None
    gancho_pendente: bool = False
    gancho_90: object = None  # Geometria esquemática; não é aprovação normativa.


@dataclass(frozen=True)
class GeometriaGancho90:
    """Eixo esquemático de uma dobra de 90 graus, em centímetros.

    O comprimento desenvolvido considera a tangente horizontal, o arco
    circular e a ponta vertical. Não substitui o detalhamento do nó.
    """
    x_inicio_cm: float
    x_tangencia_cm: float
    x_vertical_cm: float
    raio_eixo_cm: float
    ponta_vertical_cm: float
    comprimento_desenvolvido_cm: float
    sentido_vertical: int


def _barra_referencia(modelo, x_face):
    """Seleciona uma barra de seção cheia próxima à face do apoio."""
    candidatas = []
    for barra in modelo.barras:
        if getattr(barra, 'ficticia', False):
            continue
        ni, nj = barra.noi, barra.noj
        if abs(float(ni.y)) > 1e-6 or abs(float(nj.y)) > 1e-6:
            continue
        a, b = sorted((float(ni.x), float(nj.x)))
        if b - a < 1e-8:
            continue
        distancia = min(abs(a - x_face), abs(b - x_face))
        candidatas.append((distancia, barra))
    if not candidatas:
        raise ValueError('Nenhuma barra de seção cheia para avaliar ancoragem nos apoios.')
    return min(candidatas, key=lambda item: item[0])[1]


def ancorar_nos_apoios(modelo, armaduras, boa_aderencia=True, tolerancia_cm=1e-4,
                      propor_gancho=True):
    """Prolonga barras que alcançam a face do apoio, sem exceder sua face externa.

    A área calculada, quando disponível, é usada na razão As,calc/As,ef.
    Para faixas escalonadas sem essa informação, usa razão 1,0 como
    estimativa conservadora. Registra insuficiência sem criar ganchos fictícios.

    Returns:
        (armaduras_atualizadas, lista_de_diagnosticos)
    """
    # modelo.l é o vão LIVRE; x_esq e x_dir são as faces internas.
    # Os nós de apoio estão no centro dos pilares, não na face de ancoragem.
    x_esq = float(modelo.ap1) / 2.0
    x_dir = x_esq + float(modelo.l)
    exterior_esq = x_esq - float(modelo.ap1)
    exterior_dir = x_dir + float(modelo.ap2)
    if (not all(isfinite(v) for v in (x_esq, x_dir, exterior_esq, exterior_dir))
            or x_dir <= x_esq or modelo.ap1 <= 0 or modelo.ap2 <= 0):
        raise ValueError('Geometria dos apoios inválida.')
    # Reserva geométrica na face externa do pilar (modelo.cob).
    cobrimento = float(getattr(modelo, 'cob', 0.0))
    if not isfinite(cobrimento) or cobrimento < 0:
        raise ValueError('Cobrimento longitudinal inválido.')
    referencias = {
        'esquerdo': _barra_referencia(modelo, x_esq),
        'direito': _barra_referencia(modelo, x_dir),
    }
    atualizadas, diagnosticos = [], []
    for p in armaduras:
        inicio, fim = float(p.inicio_cm), float(p.fim_cm)
        if not inicio < fim:
            raise ValueError(f'Faixa inválida: {p.posicao}')
        as_ef = float(p.as_fornecida_cm2)
        as_calc = float(getattr(p, 'as_necessaria_cm2', as_ef))
        if as_calc <= 0:
            as_calc = as_ef  # Porta-estribo: avaliação conservadora provisória.
        as_calc = min(as_calc, as_ef)
        for lado, face, externo in (
            ('esquerdo', x_esq, exterior_esq),
            ('direito', x_dir, exterior_dir),
        ):
            if lado == 'esquerdo':
                toca = inicio <= face + tolerancia_cm and fim > face
                disponivel = max(0.0, face - externo)
            else:
                toca = fim >= face - tolerancia_cm and inicio < face
                disponivel = max(0.0, externo - face)
            if not toca:
                continue
            lb = comprimento_ancoragem_necessario(
                referencias[lado], p.bitola_mm, as_calc, as_ef,
                boa_aderencia=boa_aderencia)
            # O apoio começa na face da viga; não supomos concreto além da face externa.
            # Reserva geométrica preliminar de cobrimento na face externa.
            # O valor real depende da posição da barra e do detalhamento.
            disponivel = max(0.0, disponivel - cobrimento)
            # Na insuficiência, ocupar todo o espaço disponível, sem violar
            # o cobrimento. Não transformar o aviso em aprovação normativa.
            # Se a barra reta não couber, propor redução normativa de
            # comprimento para barra tracionada com gancho (alpha=0,7).
            # IMPORTANTE: não comprova viabilidade da dobra no pilar.
            lb_gancho = None
            geometria_gancho = None
            tipo = 'reta'
            pendente = False
            if propor_gancho and disponivel + tolerancia_cm < lb:
                lb_gancho = comprimento_ancoragem_necessario(
                    referencias[lado], p.bitola_mm, as_calc, as_ef,
                    boa_aderencia=boa_aderencia, alpha=0.7)
                tipo = 'gancho proposto'
                pendente = True
            lb_alvo = lb_gancho if lb_gancho is not None else lb
            adotado = min(lb_alvo, disponivel)
            if pendente:
                # A projeção horizontal é limitada ao comprimento requerido
                # para a alternativa com gancho. Quando insuficiente, ocupa
                # o espaço disponível até o cobrimento da face externa.
                # O coeficiente alpha não comprova a geometria do gancho.
                phi_cm = float(p.bitola_mm) / 10.0
                raio = 3.0 * phi_cm  # provisório: eixo da barra
                ponta = 10.0 * phi_cm  # provisório: ponta vertical
                projecao = min(lb_gancho, disponivel)
                if projecao >= raio:
                    x_tangencia = (face - (projecao - raio)
                                   if lado == 'esquerdo' else
                                   face + (projecao - raio))
                    x_vertical = (face - projecao if lado == 'esquerdo'
                                  else face + projecao)
                    sentido = 1 if p.face == 'inferior' else -1
                    trecho = projecao - raio
                    desenvolvido = trecho + pi * raio / 2 + ponta
                    geometria_gancho = GeometriaGancho90(
                        face, x_tangencia, x_vertical, raio, ponta,
                        desenvolvido, sentido)
                    adotado = projecao
                else:
                    tipo = 'gancho inviável geometricamente'
                    pendente = False
                    adotado = disponivel
            if lado == 'esquerdo':
                inicio = min(inicio, face - adotado)
            else:
                fim = max(fim, face + adotado)
            # O gancho proposto NÃO recebe status de ancoragem atendida:
            # faltam verificações de geometria, dobra e condições normativas.
            atende = (disponivel + tolerancia_cm >= lb)
            diagnosticos.append(DiagnosticoAncoragemApoio(
                p.posicao, lado, p.face, lb, disponivel, adotado, atende,
                ('Espaço geométrico suficiente para barra reta; '
                 'verificar regras específicas do apoio.' if atende else
                 'GANCHO PROPOSTO: verificar diâmetro de dobra, perna, cobrimento, confinamento e regras do apoio.'
                 if pendente else 'INSUFICIENTE: estudar solução de ancoragem específica.'),
                tipo, lb_gancho, pendente, geometria_gancho
            ))
        mudancas = {'inicio_cm': inicio, 'fim_cm': fim}
        if hasattr(p, 'comprimento_geometrico_cm'):
            mudancas['comprimento_geometrico_cm'] = fim - inicio
        atualizadas.append(replace(p, **mudancas))
    return atualizadas, diagnosticos


def relatorio_ancoragem_apoios(diagnosticos):
    """Emite no terminal as verificações efetuadas e suas insuficiências."""
    if not diagnosticos:
        print('[ANCORAGEM APOIOS] Nenhuma barra alcançou as faces dos apoios; '
              'não houve verificação. Confira a continuidade das barras.')
        return
    print(f'[ANCORAGEM APOIOS] {len(diagnosticos)} extremidade(s) verificada(s):')
    for d in diagnosticos:
        estado = ('OK GEOMÉTRICO' if d.atendida_geometricamente else
                  'GANCHO PROPOSTO (NÃO VERIFICADO)' if d.gancho_pendente
                  else 'INSUFICIENTE')
        deficit = max(0.0, d.lb_necessario_cm - d.comprimento_disponivel_cm)
        print(f'  {estado}: {d.posicao} ({d.face}), apoio {d.lado}: '
              f'lb,nec={d.lb_necessario_cm:.1f} cm; '
              f'disponível={d.comprimento_disponivel_cm:.1f} cm; '
              f'déficit reta={deficit:.1f} cm' +
              (f'; lb,nec gancho={d.lb_gancho_cm:.1f} cm; '
               f'gancho 90 esquemático={d.gancho_90 is not None}; '
               f'déficit gancho={max(0.0, d.lb_gancho_cm - d.comprimento_disponivel_cm):.1f} cm'
               if d.lb_gancho_cm is not None else ''))
    print('  Nota: verificação geométrica preliminar, não substitui '
          'as exigências normativas de ancoragem em apoios.')
