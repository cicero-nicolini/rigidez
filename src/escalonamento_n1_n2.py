"""Escalonamento TEÓRICO independente das armaduras N1 e N2.

Calcula opcionalmente a ancoragem necessária de barras retas, mas não
constitui detalhamento executivo normativo completo. Permite decalagem geométrica parametrizada da demanda. Mantém a bitola e a quantidade total propostas pela rotina original.
Unidades: cm, cm². Resultados majorados devem estar selecionados no modelo.
"""
from dataclasses import dataclass
from math import ceil, isfinite, pi
from .armaduras_longitudinais import propor_longitudinais


@dataclass(frozen=True)
class FaixaEscalonada:
    posicao: str
    face: str
    quantidade: int
    bitola_mm: float
    inicio_cm: float
    fim_cm: float
    nivel: int  # 0: barras de base; 1..N: grupos adicionais

    @property
    def as_fornecida_cm2(self):
        return self.quantidade * pi * (self.bitola_mm / 10)**2 / 4


def _amostras(modelo, x0, x1):
    """Áreas de aço nas extremidades de elementos horizontais de seção cheia."""
    segmentos = []
    for b in modelo.barras:
        if getattr(b, 'ficticia', False):
            continue
        a, z = b.noi, b.noj
        xa, xb = float(a.x), float(z.x)
        if (abs(float(a.y)) > 1e-6 or abs(float(z.y)) > 1e-6
                or abs(xb-xa) < 1e-8 or min(xa, xb) < x0-1e-6
                or max(xa, xb) > x1+1e-6):
            continue
        fl = getattr(modelo, "resultados_calculo", {}).get("fl", {}).get(b.num, b.fl)
        if len(fl) < 6:
            raise ValueError('Vetor de esforços inválido.')
        valores = []
        for k, m in enumerate((-float(fl[2]), float(fl[5]))):
            as1, as2 = float(b.ASL[k]), float(b.ASL2[k])
            if not all(isfinite(v) and v >= -1e-9 for v in (as1, as2)) or not isfinite(m):
                raise ValueError('Esforços ou armaduras inválidos.')
            valores.append((max(0, as1), max(0, as2)) if m >= 0
                          else (max(0, as2), max(0, as1)))
        if xa > xb:
            xa, xb = xb, xa
            valores.reverse()
        segmentos.append((xa, xb, valores[0], valores[1]))
    return sorted(segmentos, key=lambda v: (v[0], v[1]))


def _intervalos_acima(segmentos, face_idx, limiar, x0, x1):
    """Intervalos onde As(x) > limiar, sem atravessar lacunas da malha."""
    intervalos = []
    for xa, xb, va, vb in segmentos:
        ya, yb = va[face_idx], vb[face_idx]
        if max(ya, yb) <= limiar + 1e-9:
            continue
        if min(ya, yb) > limiar + 1e-9:
            esq, dir = xa, xb
        else:
            corte = xa + (limiar - ya) * (xb - xa) / (yb - ya)
            esq, dir = (corte, xb) if ya <= limiar else (xa, corte)
        esq, dir = max(x0, esq), min(x1, dir)
        if dir - esq <= 1e-8:
            continue
        if intervalos and esq <= intervalos[-1][1] + 1e-6:
            intervalos[-1] = (intervalos[-1][0], max(dir, intervalos[-1][1]))
        else:
            intervalos.append((esq, dir))
    return intervalos


def _decalar_intervalos(intervalos, a, b, decalagem_cm):
    """Amplia e une faixas de demanda dentro de uma região de seção cheia.

    Args:
        intervalos: Intervalos teóricos (início, fim) em centímetros.
        a: Limite esquerdo da região, em centímetros.
        b: Limite direito da região, em centímetros.
        decalagem_cm: Distância de decalagem horizontal em centímetros.

    Returns:
        Lista ordenada de intervalos decalados, sem sobreposição.

    Notes:
        Não incorpora comprimento de ancoragem. A ampliação é bilateral,
        representando a envoltória de demanda deslocada em ±a_l.
    """
    resultado = []
    for inicio, fim in sorted(intervalos):
        esquerda = max(a, inicio - decalagem_cm)
        direita = min(b, fim + decalagem_cm)
        if direita - esquerda <= 1e-8:
            continue
        if resultado and esquerda <= resultado[-1][1] + 1e-6:
            resultado[-1] = (resultado[-1][0], max(direita, resultado[-1][1]))
        else:
            resultado.append((esquerda, direita))
    return resultado


def _subtrair(intervalo, ocupados):
    """Partes de intervalo sem demanda resistente."""
    partes = [intervalo]
    for a, b in ocupados:
        novas = []
        for x, y in partes:
            if a > x + 1e-8:
                novas.append((x, min(a, y)))
            if b < y - 1e-8:
                novas.append((max(b, x), y))
        partes = [(x, y) for x, y in novas if y-x > 1e-8]
    return partes


def distancia_decalagem(barra, automatica=True, manual_cm=None):
    """Distância geométrica a_l (cm) compartilhada por diagrama e armaduras.

    Usa d, modelo de cortante e inclinações armazenadas na barra.
    Esta aproximação geométrica ainda requer conferência normativa.
    """
    from math import radians, tan

    if manual_cm is not None:
        valor = float(manual_cm)
        if not isfinite(valor) or valor < 0:
            raise ValueError('decalagem_cm deve ser finita e não negativa.')
        return valor
    if not automatica:
        return 0.0
    d = float(getattr(barra.secao, 'd', 0.0))
    theta = 45.0 if getattr(barra, 'modelo', 1) == 1 else float(barra.θ)
    alfa = float(getattr(barra, 'α', 90.0))
    if not all(isfinite(v) for v in (d, theta, alfa)):
        return 0.0
    if d <= 0 or not (0 < theta < 90 and 0 < alfa <= 90):
        return 0.0
    return max(0.0, d * (1 / tan(radians(theta)) -
                         1 / tan(radians(alfa))) / 2)


def _decalagem_trecho(modelo, a, b, automatica, manual_cm):
    """Maior a_l das barras de seção cheia pertencentes ao trecho."""
    distancias = []
    for barra in modelo.barras:
        if getattr(barra, 'ficticia', False):
            continue
        noi, noj = barra.noi, barra.noj
        xa, xb = sorted((float(noi.x), float(noj.x)))
        if (abs(float(noi.y)) > 1e-6 or abs(float(noj.y)) > 1e-6
                or xb - xa < 1e-8 or xa < a - 1e-6 or xb > b + 1e-6):
            continue
        distancias.append(distancia_decalagem(barra, automatica, manual_cm))
    return max(distancias, default=0.0)


def comprimento_ancoragem_basico(barra, bitola_mm, boa_aderencia=True):
    """Estima l_b de barra reta tracionada, em cm (NBR 6118).

    Usa f_bd = eta1 * eta2 * eta3 * f_ctd e l_b = phi/4 * f_yd/f_bd.
    Requer que as propriedades de aço e concreto estejam em unidades coerentes.
    Não inclui ganchos, emendas ou confinamento.

    Args:
        barra: Elemento com concreto.fctd e aco.fyd.
        bitola_mm: Diâmetro nominal da barra longitudinal, em mm.
        boa_aderencia: True para região de boa aderência; False para má.

    Returns:
        Comprimento básico de ancoragem em cm.

    Raises:
        ValueError: Se faltarem propriedades ou valores válidos.
    """
    phi = float(bitola_mm) / 10.0
    if phi <= 0 or not isfinite(phi):
        raise ValueError('Bitola de ancoragem inválida.')
    concreto, aco = getattr(barra, 'concreto', None), getattr(barra, 'aco', None)
    fctd = float(getattr(concreto, 'fctd', 0.0))
    fyd = float(getattr(aco, 'fyd', 0.0))
    if not all(isfinite(v) and v > 0 for v in (fctd, fyd)):
        raise ValueError('Ancoragem requer concreto.fctd e aco.fyd positivos.')
    eta1 = 2.25  # Barras nervuradas; confirmar a classe real da barra.
    eta2 = 1.0 if boa_aderencia else 0.7
    eta3 = 1.0 if bitola_mm <= 32 else (132.0 - bitola_mm) / 100.0
    if eta3 <= 0:
        raise ValueError('Bitola fora da faixa suportada para ancoragem.')
    fbd = eta1 * eta2 * eta3 * fctd
    return phi * fyd / (4.0 * fbd)


def comprimento_ancoragem_necessario(barra, bitola_mm, as_calc_cm2,
                                    as_ef_cm2, boa_aderencia=True,
                                    alpha=1.0):
    """Calcula l_b,nec de barras retas tracionadas, em cm.

    l_b,nec = max(alpha * l_b * As,calc / As,ef; l_b,min), com
    l_b,min = max(0,3*l_b; 10*phi; 10 cm).

    Args:
        barra: Elemento com materiais e propriedades de resistência.
        bitola_mm: Diâmetro da barra, em milímetros.
        as_calc_cm2: Área necessária atribuída ao grupo ancorado, em cm².
        as_ef_cm2: Área efetiva desse mesmo grupo, em cm².
        boa_aderencia: Condição de aderência definida pelo projetista.
        alpha: Coeficiente de ancoragem (1,0 para barra reta tracionada).

    Returns:
        Comprimento necessário de ancoragem, em centímetros.

    Raises:
        ValueError: Se os parâmetros forem inválidos.
    """
    as_calc = float(as_calc_cm2)
    as_ef = float(as_ef_cm2)
    alpha = float(alpha)
    if not all(isfinite(v) for v in (as_calc, as_ef, alpha)):
        raise ValueError('Áreas e coeficiente de ancoragem devem ser finitos.')
    if as_calc < 0 or as_ef <= 0 or alpha <= 0:
        raise ValueError('Áreas ou coeficiente de ancoragem inválidos.')
    if as_calc > as_ef + 1e-6:
        raise ValueError('As,calc não pode superar As,ef do grupo ancorado.')

    lb = comprimento_ancoragem_basico(barra, bitola_mm, boa_aderencia)
    phi_cm = float(bitola_mm) / 10.0
    lb_min = max(0.3 * lb, 10.0 * phi_cm, 10.0)
    return max(alpha * lb * as_calc / as_ef, lb_min)


def _barras_trecho(modelo, a, b):
    """Seleciona elementos horizontais de seção cheia de um trecho."""
    barras = []
    for barra in modelo.barras:
        if getattr(barra, 'ficticia', False):
            continue
        noi, noj = barra.noi, barra.noj
        xa, xb = sorted((float(noi.x), float(noj.x)))
        if (abs(float(noi.y)) > 1e-6 or abs(float(noj.y)) > 1e-6
                or xb - xa < 1e-8 or xa < a - 1e-6 or xb > b + 1e-6):
            continue
        barras.append(barra)
    return barras


def _ancoragem_necessaria_trecho(barras, bitola_mm, as_calc, as_ef,
                                 boa_aderencia):
    """Adota o maior l_b,nec entre os elementos da região analisada."""
    if not barras:
        raise ValueError('Sem elementos para calcular ancoragem do trecho.')
    return max(comprimento_ancoragem_necessario(
        barra, bitola_mm, as_calc, as_ef, boa_aderencia)
        for barra in barras)


def _demanda_maxima(segmentos, face_idx, inicio, fim):
    """Maior As requerida no intervalo, com interpolação nas extremidades."""
    valores = []
    for xa, xb, va, vb in segmentos:
        esq, dir = max(inicio, xa), min(fim, xb)
        if dir < esq - 1e-8:
            continue
        for x in (esq, dir):
            t = (x - xa) / (xb - xa)
            valores.append(va[face_idx] + t * (vb[face_idx] - va[face_idx]))
    return max(valores, default=0.0)


def _intervalos_ancorados(intervalos, a, b, segmentos, face_idx,
                          as_previa, as_grupo, barras, bitola_mm,
                          boa_aderencia, habilitada):
    """Decalação já aplicada; acrescenta l_b,nec a cada faixa resistente.

    A área calculada atribuída ao grupo é a demanda máxima da faixa
    menos a capacidade dos grupos anteriores. A área efetiva corresponde
    somente às barras deste grupo. Essa aproximação não substitui a
    verificação local de tensões nas seções de interrupção.
    """
    resultado = []
    for inicio, fim in intervalos:
        if habilitada:
            demanda = _demanda_maxima(segmentos, face_idx, inicio, fim)
            as_calc = min(as_grupo, max(0.0, demanda - as_previa))
            lb_nec = _ancoragem_necessaria_trecho(
                barras, bitola_mm, as_calc, as_grupo, boa_aderencia)
        else:
            lb_nec = 0.0
        esquerda = max(a, inicio - lb_nec)
        direita = min(b, fim + lb_nec)
        if resultado and esquerda <= resultado[-1][1] + 1e-6:
            resultado[-1] = (resultado[-1][0], max(direita, resultado[-1][1]))
        else:
            resultado.append((esquerda, direita))
    return resultado


def propor_escalonamento(modelo, niveis_n1=2, niveis_n2=2,
                          barras_base_n1=2, barras_base_n2=2, decalagem_cm=None,
                          decalagem_automatica=True, ancoragem_automatica=False,
                          boa_aderencia=True):
    """Escalonamento teórico por face em trechos cheios.

    Onde As requerida é zero, propõe 2 porta-estribos com bitola igual
    à dos estribos da seção cheia correspondente. Não verifica ancoragem,
    continuidade física nem comprimentos executivos.

    Args:
        decalagem_cm: Ampliação bilateral da envoltória de demanda, em cm.
        ancoragem_automatica: Acrescenta l_b,nec a cada extremidade livre das faixas
            resistentes, após a decalagem, sem ultrapassar o trecho.
        boa_aderencia: Condição de aderência das barras (deve ser verificada
            pelo projetista).
    """
    from .armaduras_transversais import propor_transversais

    for nome, valor in (('niveis_n1', niveis_n1), ('niveis_n2', niveis_n2)):
        if type(valor) is not int or valor < 0:
            raise ValueError(f'{nome} deve ser inteiro não negativo.')
    if decalagem_cm is not None:
        decalagem_cm = distancia_decalagem(None, manual_cm=decalagem_cm)
    originais = propor_longitudinais(modelo)
    estribos = propor_transversais(modelo)
    x0, x1 = float(modelo.ap1)/2, float(modelo.ap1)/2+float(modelo.l)
    segmentos = _amostras(modelo, x0, x1)
    if not segmentos:
        raise ValueError('Nenhum elemento de seção cheia encontrado.')
    saida = []
    for original in originais:
        limite = niveis_n1 if original.posicao == 'N1' else niveis_n2
        base_min = barras_base_n1 if original.posicao == 'N1' else barras_base_n2
        face_idx = 0 if original.face == 'inferior' else 1
        n = original.quantidade
        area_barra = pi * (original.bitola_mm/10)**2/4
        for trecho in estribos:
            a, b = trecho.inicio_cm, trecho.fim_cm
            decalagem_trecho = _decalagem_trecho(
                modelo, a, b, decalagem_automatica, decalagem_cm)
            locais = [seg for seg in segmentos if seg[0] >= a-1e-6 and seg[1] <= b+1e-6]
            barras_trecho = (_barras_trecho(modelo, a, b)
                             if ancoragem_automatica else [])
            demanda_decalada = _decalar_intervalos(
                _intervalos_acima(locais, face_idx, 1e-9, a, b),
                a, b, decalagem_trecho)
            positivos = _intervalos_ancorados(
                demanda_decalada, a, b, locais, face_idx,
                0.0, n * area_barra, barras_trecho,
                original.bitola_mm, boa_aderencia, ancoragem_automatica)
            if not positivos:
                for ini, fim in _subtrair((a, b), positivos):
                    saida.append(FaixaEscalonada(original.posicao, original.face,
                                2, trecho.bitola_mm, ini, fim, -1))
                continue
            if limite == 0 or n <= base_min:
                for ini, fim in _subtrair((a, b), positivos):
                    saida.append(FaixaEscalonada(original.posicao, original.face,
                                2, trecho.bitola_mm, ini, fim, -1))
                # No modo 0, mantém a proposta resistente original, sem
                # escalonamento, dentro de cada trecho de seção cheia.
                for ini, fim in positivos:
                    saida.append(FaixaEscalonada(original.posicao, original.face,
                                n, original.bitola_mm, ini, fim, 0))
                continue
            base = min(n, max(2, int(base_min)))
            positivos = _intervalos_ancorados(
                demanda_decalada, a, b, locais, face_idx,
                0.0, base * area_barra, barras_trecho,
                original.bitola_mm, boa_aderencia, ancoragem_automatica)
            for ini, fim in _subtrair((a, b), positivos):
                saida.append(FaixaEscalonada(original.posicao, original.face,
                            2, trecho.bitola_mm, ini, fim, -1))
            for ini, fim in positivos:
                saida.append(FaixaEscalonada(original.posicao, original.face,
                            base, original.bitola_mm, ini, fim, 0))
            excedentes = n-base
            grupos = min(limite, excedentes)
            tamanhos = [excedentes//grupos + (i < excedentes%grupos)
                        for i in range(grupos)]
            acumulado = base
            for nivel, tamanho in enumerate(tamanhos, 1):
                demanda_nivel = _decalar_intervalos(
                    _intervalos_acima(locais, face_idx,
                                      acumulado * area_barra, a, b),
                    a, b, decalagem_trecho)
                faixas_nivel = _intervalos_ancorados(
                    demanda_nivel, a, b, locais, face_idx,
                    acumulado * area_barra, tamanho * area_barra,
                    barras_trecho, original.bitola_mm,
                    boa_aderencia, ancoragem_automatica)
                for ini, fim in faixas_nivel:
                    saida.append(FaixaEscalonada(original.posicao, original.face,
                                 tamanho, original.bitola_mm, ini, fim, nivel))
                acumulado += tamanho
    return saida


def continuidade_nos_furos(modelo, faixas, tolerancia=1e-5):
    """Estende faixas N1/N2 incidentes a uma abertura através do banzo externo.

    Retorna (faixas_continuas, areas_por_furo). Não calcula ancoragens.
    Uma faixa que chega a qualquer borda do furo atravessa toda a abertura;
    faixas coincidentes dos dois lados são unidas sem duplicar aço.
    """
    from dataclasses import replace
    from .armaduras_banzos import _zonas
    zonas = list(_zonas(modelo))
    resultado = list(faixas)
    for _, xa, xb, _, _ in zonas:
        novos = []
        for f in resultado:
            # A faixa deve alcançar a borda do furo; não se prolongam barras
            # que terminam afastadas dele (isso exigiria verificação de ancoragem).
            esquerda = f.inicio_cm < xa-tolerancia and f.fim_cm >= xa-tolerancia
            direita = f.fim_cm > xb+tolerancia and f.inicio_cm <= xb+tolerancia
            if esquerda or direita:
                f = replace(f, inicio_cm=min(f.inicio_cm, xa),
                            fim_cm=max(f.fim_cm, xb))
            novos.append(f)
        resultado = novos
    # Une trechos idênticos fisicamente, inclusive quando coincidem no furo.
    agrupados = {}
    for f in resultado:
        chave = (f.posicao, f.face, f.quantidade, f.bitola_mm, f.nivel)
        agrupados.setdefault(chave, []).append(f)
    unidos = []
    for grupo in agrupados.values():
        grupo.sort(key=lambda f: (f.inicio_cm, f.fim_cm))
        atual = grupo[0]
        for f in grupo[1:]:
            if f.inicio_cm <= atual.fim_cm + tolerancia:
                atual = replace(atual, fim_cm=max(atual.fim_cm, f.fim_cm))
            else:
                unidos.append(atual)
                atual = f
        unidos.append(atual)
    # Nas aberturas, evitar porta-estribos duplicados na mesma face.
    # A escolha considera a bitola do estribo do banzo correspondente.
    from .armaduras_estribos_banzos import propor_estribos_banzos
    estribos_banzos = propor_estribos_banzos(modelo)
    bitolas_minimas = {
        (e.indice_furo, 'N1' if e.banzo == 'inferior' else 'N2'):
        float(e.bitola_mm) for e in estribos_banzos
    }
    for indice, xa, xb, _, _ in zonas:
        for posicao in ('N1', 'N2'):
            candidatos = [f for f in unidos
                          if f.posicao == posicao and f.nivel == -1
                          and f.inicio_cm <= xa + tolerancia
                          and f.fim_cm >= xb - tolerancia]
            if len(candidatos) < 2:
                continue
            minimo = bitolas_minimas.get((indice, posicao), 0.0)
            elegiveis = [f for f in candidatos
                         if f.bitola_mm + tolerancia >= minimo]
            if not elegiveis:
                # Sem bitola adequada, não eliminar barras sem substituição.
                continue
            escolhido = max(elegiveis,
                            key=lambda f: (f.bitola_mm, f.quantidade))
            novos = []
            for f in unidos:
                if f not in candidatos or f is escolhido:
                    novos.append(f)
                    continue
                # Preserva as partes necessárias fora do furo.
                if f.inicio_cm < xa - tolerancia:
                    novos.append(replace(f, fim_cm=xa))
                if f.fim_cm > xb + tolerancia:
                    novos.append(replace(f, inicio_cm=xb))
            unidos = novos

    areas = {}
    for indice, xa, xb, _, _ in zonas:
        areas[indice] = {}
        for posicao in ('N1', 'N2'):
            # Só contam faixas que cobrem integralmente a abertura.
            areas[indice][posicao] = sum(
                f.as_fornecida_cm2 for f in unidos
                if f.posicao == posicao and f.inicio_cm <= xa+tolerancia
                and f.fim_cm >= xb-tolerancia)
    return unidos, areas
