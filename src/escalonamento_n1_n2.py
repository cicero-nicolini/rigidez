"""Escalonamento TEÓRICO independente das armaduras N1 e N2.

Não calcula ancoragens, decalagem, comprimentos de corte ou detalhamento
normativo. Mantém a bitola e a quantidade total propostas pela rotina original.
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


def propor_escalonamento(modelo, niveis_n1=2, niveis_n2=2,
                          barras_base_n1=2, barras_base_n2=2):
    """Escalonamento teórico por face em trechos cheios.

    Onde As requerida é zero, propõe 2 porta-estribos com bitola igual
    à dos estribos da seção cheia correspondente. Não verifica ancoragem,
    continuidade física nem comprimentos executivos.
    """
    from .armaduras_transversais import propor_transversais

    for nome, valor in (('niveis_n1', niveis_n1), ('niveis_n2', niveis_n2)):
        if type(valor) is not int or valor < 0:
            raise ValueError(f'{nome} deve ser inteiro não negativo.')
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
            locais = [seg for seg in segmentos if seg[0] >= a-1e-6 and seg[1] <= b+1e-6]
            positivos = _intervalos_acima(locais, face_idx, 1e-9, a, b)
            # Porta-estribos são propostos apenas onde As da face é nula.
            for ini, fim in _subtrair((a, b), positivos):
                saida.append(FaixaEscalonada(original.posicao, original.face,
                           2, trecho.bitola_mm, ini, fim, -1))
            if not positivos:
                continue
            if limite == 0 or n <= base_min:
                # No modo 0, mantém a proposta resistente original, sem
                # escalonamento, dentro de cada trecho de seção cheia.
                for ini, fim in positivos:
                    saida.append(FaixaEscalonada(original.posicao, original.face,
                                n, original.bitola_mm, ini, fim, 0))
                continue
            base = min(n, max(2, int(base_min)))
            for ini, fim in positivos:
                saida.append(FaixaEscalonada(original.posicao, original.face,
                            base, original.bitola_mm, ini, fim, 0))
            excedentes = n-base
            grupos = min(limite, excedentes)
            tamanhos = [excedentes//grupos + (i < excedentes%grupos)
                        for i in range(grupos)]
            acumulado = base
            for nivel, tamanho in enumerate(tamanhos, 1):
                for ini, fim in _intervalos_acima(locais, face_idx,
                                                   acumulado*area_barra, a, b):
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
