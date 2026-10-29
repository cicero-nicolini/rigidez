"""Visualização de modelos de pórtico plano com Matplotlib."""

# Funções legadas preservadas para compatibilidade com os cálculos.
def print_matriz(m) -> None:
    l, c = m.shape
    tamanho_colunas = []
    for i in range(c):
        maior_l = 0
        for j in range(l):
            s = f"{m[j,i]:.1f}"
            maior_l = max(len(s), maior_l)
        tamanho_colunas.append(maior_l)
    cabecalho = '\x1b[90m  '
    for i in range(c):
        cabecalho += ' ' + str(i).center(tamanho_colunas[i])
    cabecalho += '\x1b[0m'
    print(cabecalho)
    for i in range(l):
        linha = '\x1b[90m' + str(i).rjust(2) + '\x1b[0m'
        for j in range(c):
            s = f"{m[j,i]:.1f}".rjust(tamanho_colunas[j])
            linha += ' ' + s
        print(linha)


def print_vetor(v) -> None:
    def fmt(n):
        if abs(n) < 0.06 or abs(n) > 1e6:
            return f'{n:#.2g}'
        return f'{n:.1f}'
    l = v.shape[0]
    maior_l = max((len(fmt(v[i])) for i in range(l)), default=0)
    for i in range(l):
        print('\x1b[90m' + str(i).rjust(2) + '\x1b[0m' +
              ' ' + fmt(v[i]).rjust(maior_l))


class Visualizacao:
    """Exibe geometria, nós, barras e apoios com controles interativos."""

    def __init__(self, modelo):
        self.modelo = modelo
        self._controles = []  # Mantém referências aos widgets Matplotlib.

    @staticmethod
    def _trechos_geometricos(barras):
        """Agrupa barras reais consecutivas, colineares e de mesma seção."""
        from collections import defaultdict
        import numpy as np

        reais = [b for b in barras if not getattr(b, 'ficticia', False)]
        por_no = defaultdict(list)
        for b in reais:
            por_no[id(b.noi)].append(b)
            por_no[id(b.noj)].append(b)

        def compativeis(a, b, no):
            if a.secao is not b.secao:
                return False
            def outro(barra):
                return barra.noj if barra.noi is no else barra.noi
            pa, pb = outro(a), outro(b)
            va = np.array([pa.x-no.x, pa.y-no.y], dtype=float)
            vb = np.array([pb.x-no.x, pb.y-no.y], dtype=float)
            norma = np.linalg.norm(va)*np.linalg.norm(vb)
            return norma > 0 and abs((va[0] * vb[1] - va[1] * vb[0])) <= 1e-8*norma and np.dot(va, vb) < 0

        vizinhos = defaultdict(list)
        for no_id, ligadas in por_no.items():
            for i, a in enumerate(ligadas):
                for b in ligadas[i+1:]:
                    no = a.noi if id(a.noi) == no_id else a.noj
                    if compativeis(a, b, no):
                        vizinhos[id(a)].append(b)
                        vizinhos[id(b)].append(a)
        grupos, vistos = [], set()
        for b in reais:
            if id(b) in vistos:
                continue
            pilha, grupo = [b], []
            while pilha:
                atual = pilha.pop()
                if id(atual) in vistos:
                    continue
                vistos.add(id(atual))
                grupo.append(atual)
                pilha.extend(vizinhos[id(atual)])
            grupos.append(grupo)
        return grupos

    @staticmethod
    def _pontos_rotulo_trechos(grupos, valores_por_barra):
        """Retorna início/fim e extremos interiores de cada trecho.

        valores_por_barra: id(barra) -> (valor_i, valor_j).
        Os extremos são pesquisados nas extremidades de todos os elementos.
        """
        saida = []
        for grupo in grupos:
            dados = []
            for b in grupo:
                if id(b) not in valores_por_barra:
                    continue
                vi, vj = valores_por_barra[id(b)]
                dados.extend([(b.noi, float(vi), b), (b.noj, float(vj), b)])
            if not dados:
                continue
            # Ordenação pela projeção ao longo da direção do trecho.
            b0 = grupo[0]
            dx, dy = b0.noj.x-b0.noi.x, b0.noj.y-b0.noi.y
            norm = (dx*dx+dy*dy)**0.5
            ux, uy = dx/norm, dy/norm
            dados.sort(key=lambda item: (item[0].x-b0.noi.x)*ux + (item[0].y-b0.noi.y)*uy)
            # Para cada nó interno, conservam-se os valores dos dois lados
            # quando existe descontinuidade de esforços.
            minimo, maximo = dados[0][0], dados[-1][0]
            candidatos = [dados[0], dados[-1]]
            interiores = [d for d in dados if d[0] is not minimo and d[0] is not maximo]
            if interiores:
                maior = max(interiores, key=lambda d: d[1])
                menor = min(interiores, key=lambda d: d[1])
                if maior[1] > max(dados[0][1], dados[-1][1]) + 1e-9:
                    candidatos.append(maior)
                if menor[1] < min(dados[0][1], dados[-1][1]) - 1e-9:
                    candidatos.append(menor)
            vistos = set()
            for no, valor, barra in candidatos:
                chave = (id(no), round(valor, 7))
                if chave in vistos:
                    continue
                vistos.add(chave)
                saida.append((no, valor, barra))
        return saida

    def calcular_momento_barra(self, barra, n_pontos=50):
        """Interpola linearmente o momento a partir de barra.fl, sem cargas.

        Convencao de momentos internos adotada:
        M(0) = -fl[2] e M(L) = fl[5].
        O diagrama e uma interpolacao grafica entre os valores de extremidade,
        nao a distribuicao exata sob carregamento entre nos.
        """
        import numpy as np

        L = float(barra.L)
        if L <= 0:
            raise ValueError(f"Barra {barra.num} com comprimento invalido.")
        fl = np.asarray(barra.esforcos_exibicao, dtype=float).ravel()
        if fl.size < 6:
            raise ValueError(f"Barra {barra.num}: barra.fl deve ter 6 componentes.")

        x = np.linspace(0.0, L, n_pontos)
        Mi, Mj = -fl[2], fl[5]
        M = Mi + (Mj - Mi) * (x / L)
        return x, M

    def desenhar_geometria(self, mostrar_nos=False, mostrar_barras=False,
                       mostrar_apoios=True, ax=None, controles=True):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon, Rectangle, Circle
        from matplotlib.widgets import CheckButtons, Button, Slider

        nos, barras = self.modelo.nos, self.modelo.barras
        if not nos or not barras:
            raise ValueError('O modelo precisa conter nós e barras.')

        if ax is None:
            fig, ax = plt.subplots(figsize=(14, 7.5))
        else:
            fig = ax.figure
        if controles:
            fig.subplots_adjust(left=0.06, right=0.77, top=0.94, bottom=0.10)

        xs, ys = [n.x for n in nos], [n.y for n in nos]
        xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
        dx, dy = xmax - xmin, ymax - ymin
        escala = max(dx, dy, 1.0)
        tam_apoio = 0.018 * escala

        # Camada opcional: um único contorno contínuo para viga e pilares.
        # A abertura é um furo real do caminho, não um retângulo branco.
        geometria_viga = []
        parametros = ('ap1', 'ap2', 'l', 'lf', 'b', 'h', 'a', 'c')
        if all(hasattr(self.modelo, nome) for nome in parametros):
            from matplotlib.path import Path as MplPath
            from matplotlib.patches import PathPatch

            m = self.modelo
            h = float(m.h)
            ap1, ap2 = float(m.ap1), float(m.ap2)
            x0 = ap1 / 2.0
            x1 = x0 + float(m.l)
            y0, y1 = -h / 2.0, h / 2.0
            xc = x0 + float(m.lf)
            xa, xb = xc - float(m.b) / 2.0, xc + float(m.b) / 2.0
            ya, yb = y0 + float(m.c), y0 + float(m.c) + float(m.a)

            if h > 0 and ap1 > 0 and ap2 > 0 and x1 > x0:
                # Centros dos pilares coincidem com os nós extremos.
                apoio_esq = min(nos, key=lambda n: n.x)
                apoio_dir = max(nos, key=lambda n: n.x)
                xl = float(apoio_esq.x) - ap1 / 2.0
                xle = float(apoio_esq.x) + ap1 / 2.0
                xrd = float(apoio_dir.x) - ap2 / 2.0
                xr = float(apoio_dir.x) + ap2 / 2.0
                ye_inf, ye_sup = float(apoio_esq.y) - h, float(apoio_esq.y) + h
                yd_inf, yd_sup = float(apoio_dir.y) - h, float(apoio_dir.y) + h

                # Percorre apenas a fronteira EXTERNA da união.
                # Não há linhas internas separando a viga dos pilares.
                contorno = [
                    (xl, ye_inf), (xl, ye_sup),
                    (xle, ye_sup), (xle, y1),
                    (xrd, y1), (xrd, yd_sup),
                    (xr, yd_sup), (xr, yd_inf),
                    (xrd, yd_inf), (xrd, y0),
                    (xle, y0), (xle, ye_inf),
                ]

                def caminho_fechado(pontos):
                    vertices = pontos + [pontos[0]]
                    codigos = ([MplPath.MOVETO]
                               + [MplPath.LINETO] * (len(pontos) - 1)
                               + [MplPath.CLOSEPOLY])
                    return MplPath(vertices, codigos)

                caminho = caminho_fechado(contorno)
                if xle < xa < xb < xrd and y0 < ya < yb < y1:
                    # Sentido inverso ao contorno externo para criar o furo.
                    abertura = caminho_fechado([
                        (xa, ya), (xb, ya), (xb, yb), (xa, yb)
                    ])
                    caminho = MplPath.make_compound_path(caminho, abertura)

                bloco = PathPatch(
                    caminho, facecolor='#808080', edgecolor='#777777',
                    linewidth=0.8, alpha=0.20, zorder=0.2
                )
                ax.add_patch(bloco)
                geometria_viga.append(bloco)

        grupos = {'Nós': [], 'Números dos nós': [], 'Números das barras': [],
                  'Apoios': [], 'Cargas distribuídas': [], 'Forças nodais': [], 'Momentos nodais': []}

        for barra in barras:
            xi, yi = barra.noi.x, barra.noi.y
            xj, yj = barra.noj.x, barra.noj.y
            ax.plot([xi, xj], [yi, yj], color='#2866a6', lw=2.1, zorder=2)
            # Deslocamento perpendicular à barra em pontos de tela.
            vx, vy = xj - xi, yj - yi
            comprimento = (vx * vx + vy * vy) ** 0.5
            if comprimento:
                ox, oy = -vy / comprimento, vx / comprimento
            else:
                ox, oy = 0, 1
            anot = ax.annotate(str(barra.num), ((xi+xj)/2, (yi+yj)/2),
                               xytext=(9*ox, 9*oy), textcoords='offset points',
                               ha='center', va='center', fontsize=8,
                               color='#b14e20', zorder=5)
            grupos['Números das barras'].append(anot)

        pontos = ax.scatter(xs, ys, s=19, color='#153956', zorder=4)
        grupos['Nós'].append(pontos)
        for no in nos:
            anot = ax.annotate(str(no.num), (no.x, no.y),
                               xytext=(4, -12), textcoords='offset points',
                               fontsize=7.5, color='#333333', zorder=5)
            grupos['Números dos nós'].append(anot)

        for no in nos:
            tx, ty, rz = bool(no.Tx), bool(no.Ty), bool(no.Rz)
            if not (tx or ty or rz):
                continue
            artistas = []
            if tx and ty and rz:  # engaste
                artistas.append(Rectangle((no.x-tam_apoio/2, no.y-tam_apoio),
                                          tam_apoio, tam_apoio, facecolor='#999999',
                                          edgecolor='black', zorder=6))
            elif tx and ty:  # apoio articulado
                artistas.append(Polygon([(no.x, no.y),
                                          (no.x-tam_apoio, no.y-tam_apoio),
                                          (no.x+tam_apoio, no.y-tam_apoio)],
                                         closed=True, facecolor='#cccccc',
                                         edgecolor='black', zorder=6))
            elif ty and not tx:  # apoio móvel vertical
                artistas.append(Polygon([(no.x, no.y),
                                          (no.x-tam_apoio, no.y-tam_apoio*.8),
                                          (no.x+tam_apoio, no.y-tam_apoio*.8)],
                                         closed=True, facecolor='#dddddd',
                                         edgecolor='black', zorder=6))
                for s in (-0.5, 0.5):
                    artistas.append(Circle((no.x+s*tam_apoio,
                                            no.y-1.03*tam_apoio),
                                           radius=tam_apoio*.16,
                                           facecolor='white', edgecolor='black',
                                           zorder=6))
            elif tx and not ty:  # apoio móvel horizontal
                artistas.append(Polygon([(no.x, no.y),
                                          (no.x-tam_apoio*.8, no.y-tam_apoio),
                                          (no.x-tam_apoio*.8, no.y+tam_apoio)],
                                         closed=True, facecolor='#dddddd',
                                         edgecolor='black', zorder=6))
            elif rz:
                artistas.append(Rectangle((no.x-tam_apoio*.25,
                                            no.y-tam_apoio*.25),
                                           tam_apoio*.5, tam_apoio*.5,
                                           fill=False, edgecolor='#555555', zorder=6))
            for artista in artistas:
                ax.add_patch(artista)
            grupos['Apoios'].extend(artistas)

        # Cargas distribuídas: Py local, positivo no eixo y local da barra.

        # Cargas distribuídas nos eixos locais das barras.
        # A representação gráfica não altera os carregamentos do modelo.
        from matplotlib.patches import FancyArrowPatch, Arc

        cargas = getattr(self.modelo, 'carregamentos', [])
        segmentos = []

        cor_carga = '#b83280'
        comprimento_seta = escala * 0.045
        espacamento = escala * 0.035

        for carga in cargas:
            if not all(hasattr(carga, k) for k in
                       ('w1', 'w2', 'lw', 'a', 'barra')):
                continue

            barra = carga.barra

            xi, yi = barra.noi.x, barra.noi.y
            xj, yj = barra.noj.x, barra.noj.y

            vx = xj - xi
            vy = yj - yi
            L = (vx**2 + vy**2)**0.5

            if L <= 0 or carga.lw <= 0:
                continue

            ex, ey = vx / L, vy / L
            nx, ny = -ey, ex

            inicio = max(0.0, carga.a)
            fim = min(L, carga.a + carga.lw)

            if fim <= inicio:
                continue

            def intensidade(s):
                return carga.w1 + (
                    carga.w2 - carga.w1
                ) * (s - carga.a) / carga.lw

            # Número de setas proporcional ao comprimento carregado.
            n = max(2, int((fim - inicio) / espacamento) + 1)

            for i in range(n):
                s = inicio + (fim - inicio) * i / (n - 1)
                w = intensidade(s)

                if abs(w) < 1e-12:
                    continue

                x = xi + ex * s
                y = yi + ey * s

                sinal = 1 if w > 0 else -1
                altura = comprimento_seta * sinal

                origem = (x + nx * altura, y + ny * altura)

                seta = FancyArrowPatch(
                    origem, (x, y),
                    arrowstyle='-|>',
                    mutation_scale=9,
                    linewidth=1.0,
                    color=cor_carga,
                    zorder=3
                )

                ax.add_patch(seta)
                grupos['Cargas distribuídas'].append(seta)

            # Linha que une as origens das setas.
            w_ini = intensidade(inicio)
            w_fim = intensidade(fim)

            def ponto_externo(s, w):
                x = xi + ex * s
                y = yi + ey * s
                sinal = 1 if w > 0 else -1
                return (
                    x + nx * comprimento_seta * sinal,
                    y + ny * comprimento_seta * sinal
                )

            p1 = ponto_externo(inicio, w_ini)
            p2 = ponto_externo(fim, w_fim)

            linha, = ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                color=cor_carga,
                linewidth=1.0,
                zorder=3
            )
            grupos['Cargas distribuídas'].append(linha)

            # Dados para agrupar as etiquetas das cargas.
            if abs(carga.w1 - carga.w2) < 1e-10:
                segmentos.append({
                    'x': (xi + xj) / 2,
                    'y': (yi + yj) / 2,
                    'valor': carga.w1,
                    'nx': nx,
                    'ny': ny
                })
            else:
                xm = (p1[0] + p2[0]) / 2
                ym = (p1[1] + p2[1]) / 2
                etiqueta = ax.annotate(
                    f'q: {carga.w1:g} a {carga.w2:g}',
                    (xm, ym),
                    xytext=(0, 9),
                    textcoords='offset points',
                    ha='center',
                    fontsize=8,
                    color=cor_carga
                )
                grupos['Cargas distribuídas'].append(etiqueta)

        # Agrupa etiquetas de cargas uniformes em faixas horizontais.
        # Uma etiqueta por faixa de mesma intensidade e orientação.
        grupos_etiquetas = {}

        for seg in segmentos:
            chave = (
                round(seg['y'], 5),
                round(seg['valor'], 8),
                round(seg['nx'], 5),
                round(seg['ny'], 5)
            )
            grupos_etiquetas.setdefault(chave, []).append(seg)

        for chave, itens in grupos_etiquetas.items():
            xmedio = sum(s['x'] for s in itens) / len(itens)
            ymedio = sum(s['y'] for s in itens) / len(itens)
            valor = itens[0]['valor']

            nx = itens[0]['nx']
            ny = itens[0]['ny']
            sinal = 1 if valor > 0 else -1

            posicao = (
                xmedio + nx * sinal * comprimento_seta * 1.25,
                ymedio + ny * sinal * comprimento_seta * 1.25
            )

            etiqueta = ax.annotate(
                f'q = {valor:g}',
                posicao,
                fontsize=9,
                fontweight='bold',
                color=cor_carga,
                ha='center',
                va='center',
                zorder=5
            )

            grupos['Cargas distribuídas'].append(etiqueta)

        # Forças nodais em coordenadas globais.
        for no in nos:
            for valor, ux, uy, desloc in ((no.Fx, 1, 0, (0,12)), (no.Fy, 0, 1, (12,0))):
                if abs(valor) < 1e-12:
                    continue
                sinal = 1 if valor > 0 else -1
                tamanho = escala*.065
                origem = (no.x-ux*sinal*tamanho, no.y-uy*sinal*tamanho)
                seta = FancyArrowPatch(origem, (no.x,no.y), arrowstyle='-|>',
                                        mutation_scale=12, lw=1.5,
                                        color='#198754', zorder=7)
                ax.add_patch(seta)
                grupos['Forças nodais'].append(seta)
                anot = ax.annotate(f'{valor:g}', origem, xytext=desloc,
                                   textcoords='offset points', fontsize=8,
                                   color='#198754', ha='center')
                grupos['Forças nodais'].append(anot)
            if abs(no.Mz) > 1e-12:
                raio = escala*.027
                arco = Arc((no.x,no.y), 2*raio, 2*raio, theta1=25, theta2=315,
                           color='#b8860b', lw=1.5, zorder=7)
                ax.add_patch(arco)
                grupos['Momentos nodais'].append(arco)
                sinal = 1 if no.Mz > 0 else -1
                ang = 315 if sinal > 0 else 25
                import math
                t = math.radians(ang)
                ponta = (no.x+raio*math.cos(t), no.y+raio*math.sin(t))
                tang = (-math.sin(t)*sinal, math.cos(t)*sinal)
                seta = FancyArrowPatch((ponta[0]-tang[0]*raio*.35,
                                        ponta[1]-tang[1]*raio*.35), ponta,
                                        arrowstyle='-|>', mutation_scale=11,
                                        color='#b8860b', zorder=8)
                ax.add_patch(seta)
                grupos['Momentos nodais'].append(seta)
                anot = ax.annotate(f'M={no.Mz:g}', (no.x,no.y),
                                   xytext=(12,15), textcoords='offset points',
                                   fontsize=8, color='#b8860b')
                grupos['Momentos nodais'].append(anot)

        # ============================================================
        # DIAGRAMA DE MOMENTO FLETOR
        # ============================================================
        artistas_momento = []
        etiquetas_momento = []
        dados_momento = []
        extremos = []
        for barra in barras:
            x_local, M = self.calcular_momento_barra(barra)
            xi, yi = barra.noi.x, barra.noi.y
            dx_b = barra.noj.x - xi
            dy_b = barra.noj.y - yi
            L = np.hypot(dx_b, dy_b)
            if L <= 0:
                continue
            c, s = dx_b/L, dy_b/L
            dados_momento.append((barra, xi, yi, c, s, x_local, M))
            if not getattr(barra, 'ficticia', False):
                extremos.extend([(float(M[0]), barra.num, 'i'),
                                 (float(M[-1]), barra.num, 'j')])

        Mmax = max((np.max(np.abs(d[6])) for d in dados_momento), default=0.0)
        escala_momento = 0.06 * escala / Mmax if Mmax > 1e-12 else 1.0
        objetos_momento = []
        etiquetas_extremos = []
        pontos_rotulados = set()
        extremos_ordenados = sorted(extremos, key=lambda item: item[0])
        minimo_global = extremos_ordenados[0] if extremos_ordenados else None
        maximo_global = extremos_ordenados[-1] if extremos_ordenados else None

        for barra, xi, yi, c, s, x_local, M in dados_momento:
            x_base = xi + c*x_local
            y_base = yi + s*x_local
            linha, = ax.plot([], [], color='#16803c', lw=1.8,
                             zorder=7, visible=False)
            preenchimento = ax.fill([], [], color='#16803c',
                                     alpha=0.13, zorder=3,
                                     visible=False)[0]
            artistas_momento.extend((linha, preenchimento))
            objetos_momento.append((linha, preenchimento,
                                     x_base, y_base, c, s, M))

        trechos = self._trechos_geometricos(barras)
        momentos_por_barra = {id(b): (float(m[0]), float(m[-1]))
                              for b, _, _, _, _, _, m in dados_momento}
        rotulos_m = self._pontos_rotulo_trechos(trechos, momentos_por_barra)
        for no, valor, barra in rotulos_m:
            c = (barra.noj.x-barra.noi.x)/barra.L
            sn = (barra.noj.y-barra.noi.y)/barra.L
            xb, yb = float(no.x), float(no.y)
            anot = ax.annotate(f'{valor:.3g}', (xb, yb),
                               xytext=(5, 8), textcoords='offset points',
                               fontsize=7, color='#126a32', visible=False, zorder=10)
            etiquetas_momento.append((anot, xb, yb, c, sn, valor))

        # Extremos globais continuam disponíveis no controle independente.
        if rotulos_m:
            for nome, item in (('Mmín', min(rotulos_m, key=lambda r: r[1])),
                               ('Mmáx', max(rotulos_m, key=lambda r: r[1]))):
                no, valor, barra = item
                c = (barra.noj.x-barra.noi.x)/barra.L
                sn = (barra.noj.y-barra.noi.y)/barra.L
                xb, yb = float(no.x), float(no.y)
                anot = ax.annotate(f'{nome} = {valor:.4g}', (xb, yb),
                                   xytext=(9, -15 if valor >= 0 else 12),
                                   textcoords='offset points', fontsize=8,
                                   fontweight='bold', color='#126a32',
                                   bbox=dict(boxstyle='round,pad=0.2',
                                             facecolor='white', alpha=0.85,
                                             edgecolor='none'),
                                   visible=False, zorder=12)
                etiquetas_extremos.append((anot, xb, yb, c, sn, valor))

        def atualizar_momento(fator):
            for linha, preenchimento, xb, yb, c, s, M in objetos_momento:
                xd = xb + s*fator*M
                yd = yb - c*fator*M
                linha.set_data(xd, yd)
                preenchimento.set_xy(np.column_stack((
                    np.r_[xb, xd[::-1]], np.r_[yb, yd[::-1]])))
            for anot, xb, yb, c, s, valor in (
                    etiquetas_momento + etiquetas_extremos):
                anot.xy = (xb + s*fator*valor, yb - c*fator*valor)
            fig.canvas.draw_idle()

        atualizar_momento(escala_momento)
        grupos['Momento fletor'] = artistas_momento
        grupos['Valores dos momentos'] = [item[0] for item in etiquetas_momento]
        grupos['Extremos do momento'] = [item[0] for item in etiquetas_extremos]

        if extremos:
            maximo = max(extremos, key=lambda item: item[0])
            minimo = min(extremos, key=lambda item: item[0])
            texto_extremos = (
                f'Mmáx = {maximo[0]:.4g} (barra {maximo[1]}, {maximo[2]})\n'
                f'Mmín = {minimo[0]:.4g} (barra {minimo[1]}, {minimo[2]})'
            )
        else:
            texto_extremos = 'Momentos indisponíveis'

        # ============================================================
        # DIAGRAMAS DE ESFORCO NORMAL E CORTANTE (somente barra.fl)
        # ============================================================
        diagramas_nv = {}
        for tipo, componentes, cor in (
            ('Esforço normal', (0, 3), '#7c3aed'),
            ('Esforço cortante', (1, 4), '#e07b18'),
        ):
            dados = []
            for barra in barras:
                fl = np.asarray(barra.esforcos_exibicao, dtype=float).ravel()
                if fl.size < 6:
                    raise ValueError(f'Barra {barra.num}: fl deve ter 6 componentes.')
                xi, yi = barra.noi.x, barra.noi.y
                vx = barra.noj.x - xi
                vy = barra.noj.y - yi
                L = np.hypot(vx, vy)
                if L <= 0:
                    continue
                # Esforços internos: sinais nas extremidades opostos
                # aos esforços nodais de ação do elemento no nó inicial.
                qi, qj = -float(fl[componentes[0]]), float(fl[componentes[1]])
                dados.append((barra, xi, yi, vx/L, vy/L, L, qi, qj))

            qmax = max((max(abs(d[6]), abs(d[7])) for d in dados), default=0.0)
            fator_auto = 0.06 * escala / qmax if qmax > 1e-12 else 1.0
            objetos, etiquetas, pontos_vistos = [], [], set()
            for barra, xi, yi, c, s, L, qi, qj in dados:
                xb = np.array([xi, xi + c*L])
                yb = np.array([yi, yi + s*L])
                q = np.array([qi, qj])
                linha, = ax.plot([], [], color=cor, lw=1.8, zorder=7,
                                 visible=False)
                preenchimento = ax.fill([], [], color=cor, alpha=0.13,
                                         zorder=3, visible=False)[0]
                objetos.append((linha, preenchimento, xb, yb, c, s, q))
            valores_nv = {id(b): (qi, qj)
                          for b, _, _, _, _, _, qi, qj in dados}
            for no, valor, barra in self._pontos_rotulo_trechos(trechos, valores_nv):
                c = (barra.noj.x-barra.noi.x)/barra.L
                sn = (barra.noj.y-barra.noi.y)/barra.L
                xb, yb = float(no.x), float(no.y)
                anot = ax.annotate(f'{valor:.3g}', (xb, yb),
                                   xytext=(5, 8), textcoords='offset points',
                                   fontsize=7, color=cor, visible=False, zorder=10)
                etiquetas.append((anot, xb, yb, c, sn, valor))

            def atualizar(fator, objetos=objetos, etiquetas=etiquetas):
                for linha, preenchimento, xb, yb, c, s, q in objetos:
                    xd = xb + s*fator*q
                    yd = yb - c*fator*q
                    linha.set_data(xd, yd)
                    preenchimento.set_xy(np.column_stack((
                        np.r_[xb, xd[::-1]], np.r_[yb, yd[::-1]])))
                for anot, xb, yb, c, s, valor in etiquetas:
                    anot.xy = (xb + s*fator*valor, yb - c*fator*valor)
                fig.canvas.draw_idle()

            atualizar(fator_auto)
            grupos[tipo] = [a for linha, preenchimento, *_ in objetos
                            for a in (linha, preenchimento)]
            grupos['Valores ' + ('N' if tipo == 'Esforço normal' else 'V')] = [
                item[0] for item in etiquetas]
            diagramas_nv[tipo] = {
                'atualizar': atualizar,
                'fator': fator_auto,
                'etiquetas': etiquetas,
                'nome_valores': 'Valores ' + ('N' if tipo == 'Esforço normal' else 'V'),
            }

        # ============================================================
        # ARMADURAS: pontos i e j ligados por retas.
        # ASL acompanha o sinal do momento em cada extremidade:
        # M negativo -> acima (barra horizontal i->j); M positivo -> abaixo.
        # ASL2 usa o sentido oposto ao ASL em cada ponto.
        # Um unico controle e escala para ambas.
        # ============================================================
        diagramas_armadura = {}
        configuracoes_armadura = (
            ('ASL', 'Armadura longitudinal', 'Valores armadura longitudinal',
             '#087e8b', +1),
            ('ASL2', 'Armadura longitudinal', 'Valores armadura longitudinal',
             '#1f9cc1', -1),
            ('ASW_S', 'Armadura transversal', 'Valores armadura transversal',
             '#b23a48', +1),
        )
        dados_armadura = {}
        for atributo, titulo, nome_valores, cor, sentido in configuracoes_armadura:
            dados = []
            for barra in barras:
                if getattr(barra, "ficticia", False):
                    continue
                bruto = getattr(barra, atributo, None)
                if bruto is None:
                    continue
                valores = np.asarray(bruto, dtype=float).ravel()
                if valores.size != 2:
                    raise ValueError(
                        f'Barra {barra.num}: {atributo} deve ter dois valores (i, j).')
                if not np.all(np.isfinite(valores)):
                    continue
                # ASL2 ausente ou nulo nao gera diagrama.
                if atributo == 'ASL2' and not np.any(np.abs(valores) > 1e-12):
                    continue
                xi, yi = barra.noi.x, barra.noi.y
                vx, vy = barra.noj.x - xi, barra.noj.y - yi
                L = np.hypot(vx, vy)
                if L <= 0:
                    continue
                # Sinais do momento interno nos pontos i e j, conforme
                # a convencao ja utilizada no diagrama de momento fletor.
                if atributo in ('ASL', 'ASL2'):
                    fl = np.asarray(barra.fl, dtype=float).ravel()
                    if fl.size < 6:
                        raise ValueError(
                            f'Barra {barra.num}: fl deve ter 6 componentes.')
                    momentos = np.array([-fl[2], fl[5]], dtype=float)
                    sentidos = np.sign(momentos)
                    # Em M=0, conserva o lado do outro extremo da barra;
                    # se ambos forem zero, usa o lado positivo do diagrama.
                    if sentidos[0] == 0:
                        sentidos[0] = sentidos[1] if sentidos[1] != 0 else 1.0
                    if sentidos[1] == 0:
                        sentidos[1] = sentidos[0]
                else:
                    sentidos = np.ones(2, dtype=float)
                dados.append((barra, xi, yi, vx/L, vy/L, L, valores, sentidos))
            dados_armadura[atributo] = dados

        for titulo, atributos in (
            ('Armadura longitudinal', ('ASL', 'ASL2')),
            ('Armadura transversal', ('ASW_S',)),
        ):
            maior = max((float(np.max(np.abs(item[-2])))
                         for atributo in atributos
                         for item in dados_armadura[atributo]), default=0.0)
            fator_auto = 0.06 * escala / maior if maior > 1e-12 else 1.0
            objetos, etiquetas, pontos_vistos = [], [], set()
            nome_valores = ('Valores armadura longitudinal' if
                           titulo == 'Armadura longitudinal' else
                           'Valores armadura transversal')
            grupos[titulo] = []
            grupos[nome_valores] = []
            for atributo, grupo, _, cor, sentido in configuracoes_armadura:
                if grupo != titulo:
                    continue
                for barra, xi, yi, c, sn, L, valores, sentidos_m in dados_armadura[atributo]:
                    xb = np.array([xi, xi + c*L])
                    yb = np.array([yi, yi + sn*L])
                    linha, = ax.plot([], [], color=cor, lw=1.8, marker='o',
                                     markersize=3.2, zorder=7, visible=False)
                    preenchimento = ax.fill([], [], color=cor, alpha=0.14,
                                             zorder=3, visible=False)[0]
                    objetos.append((linha, preenchimento, xb, yb, c, sn,
                                    valores, sentido*sentidos_m))
                    grupos[titulo].extend((linha, preenchimento))
            # Etiquetas apenas nas extremidades e nos extremos interiores
            # de cada trecho, separadamente para cada tipo de armadura.
            for atributo, grupo, _, cor, sentido in configuracoes_armadura:
                if grupo != titulo:
                    continue
                valores_atributo = {
                    id(barra): (float(valores[0]), float(valores[1]))
                    for barra, _, _, _, _, _, valores, _ in dados_armadura[atributo]
                }
                for no, valor, barra in self._pontos_rotulo_trechos(trechos, valores_atributo):
                    c = (barra.noj.x-barra.noi.x)/barra.L
                    sn = (barra.noj.y-barra.noi.y)/barra.L
                    momentos = np.asarray(barra.esforcos_exibicao, dtype=float)
                    momento = -momentos[2] if no is barra.noi else momentos[5]
                    sinal_ponto = float(sentido * (np.sign(momento) or 1.0)) if atributo != 'ASW_S' else 1.0
                    xb, yb = float(no.x), float(no.y)
                    chave = (atributo, round(xb, 5), round(yb, 5), round(valor, 8))
                    if chave in pontos_vistos:
                        continue
                    pontos_vistos.add(chave)
                    anot = ax.annotate(f'{valor:.3g}', (xb, yb),
                                       xytext=(5, 8), textcoords='offset points',
                                       fontsize=7, color=cor, visible=False, zorder=10)
                    etiquetas.append((anot, xb, yb, c, sn, valor, sinal_ponto))
                    grupos[nome_valores].append(anot)

            def atualizar_armadura(fator, objetos=objetos, etiquetas=etiquetas):
                for linha, preenchimento, xb, yb, c, sn, valores, sentidos in objetos:
                    xd = xb + sn*fator*valores*sentidos
                    yd = yb - c*fator*valores*sentidos
                    linha.set_data(xd, yd)
                    preenchimento.set_xy(np.column_stack((
                        np.r_[xb, xd[::-1]], np.r_[yb, yd[::-1]])))
                for anot, x0, y0, c, sn, valor, sentido in etiquetas:
                    anot.xy = (x0 + sentido*sn*fator*valor,
                               y0 - sentido*c*fator*valor)
                fig.canvas.draw_idle()

            atualizar_armadura(fator_auto)
            diagramas_armadura[titulo] = {
                'atualizar': atualizar_armadura,
                'fator': fator_auto,
                'nome_valores': nome_valores,
            }

        # ============================================================
        # DEFORMADA DA ESTRUTURA
        # ============================================================

        linhas_deformada = []
        dados_deformada = []
        escala_auto = 1.0
        slider = None

        portico = getattr(self.modelo, "portico", None)

        if portico is not None and hasattr(portico, "u"):
            u = np.asarray(portico.u, dtype=float).ravel()

            if len(u) != 3 * len(nos):
                raise ValueError(
                    "O vetor de deslocamentos tem dimensão incompatível."
                )

            if not np.all(np.isfinite(u)):
                raise ValueError(
                    "O vetor de deslocamentos contém valores inválidos."
                )

            # Associação entre nós e posições no vetor global.
            # Pressupõe a mesma ordem utilizada na montagem da matriz.
            indice_no = {
                id(no): i for i, no in enumerate(nos)
            }

            def deslocamento_no(no):
                i = 3 * indice_no[id(no)]
                return u[i:i + 3]

            max_u = max(
                np.hypot(*deslocamento_no(no)[:2])
                for no in nos
            )

            if max_u > 1e-14:
                escala_auto = 0.05 * escala / max_u

            for barra in barras:
                xi, yi = barra.noi.x, barra.noi.y
                xj, yj = barra.noj.x, barra.noj.y

                dx_b = xj - xi
                dy_b = yj - yi
                L = np.hypot(dx_b, dy_b)

                if L <= 0:
                    continue

                c = dx_b / L
                s = dy_b / L

                ui = deslocamento_no(barra.noi)
                uj = deslocamento_no(barra.noj)

                # Deslocamentos globais -> locais.
                u1 = c * ui[0] + s * ui[1]
                v1 = -s * ui[0] + c * ui[1]
                t1 = ui[2]

                u2 = c * uj[0] + s * uj[1]
                v2 = -s * uj[0] + c * uj[1]
                t2 = uj[2]

                r = np.linspace(0, 1, 40)

                # Interpolação axial linear.
                ua = (1-r)*u1 + r*u2

                # Interpolação transversal de Hermite.
                H1 = 1 - 3*r**2 + 2*r**3
                H2 = L*(r - 2*r**2 + r**3)
                H3 = 3*r**2 - 2*r**3
                H4 = L*(-r**2 + r**3)

                va = H1*v1 + H2*t1 + H3*v2 + H4*t2

                dados_deformada.append(
                    (xi, yi, c, s, L, r, ua, va)
                )

                linha, = ax.plot(
                    [], [],
                    color="crimson",
                    linewidth=2.0,
                    zorder=8,
                    visible=False
                )

                linhas_deformada.append(linha)

        def atualizar_deformada(fator):
            for linha, dados in zip(
                linhas_deformada, dados_deformada
            ):
                xi, yi, c, s, L, r, ua, va = dados

                xl = L*r + fator*ua
                yl = fator*va

                x = xi + c*xl - s*yl
                y = yi + s*xl + c*yl

                linha.set_data(x, y)

            fig.canvas.draw_idle()

        if linhas_deformada:
            atualizar_deformada(escala_auto)



        grupos['Geometria da viga'] = geometria_viga

        # Os rótulos numéricos são controlados por uma única opção "Valores".
        # As coleções de artistas continuam separadas internamente.
        estados = {
            'Geometria da viga': False,
            'Nós': True,
            'Números dos nós': mostrar_nos,
            'Números das barras': mostrar_barras,
            'Apoios': mostrar_apoios,
            'Cargas distribuídas': True,
            'Forças nodais': True,
            'Momentos nodais': True,
            'Deformada': False,
            'Momento fletor': False,
            'Esforço normal': False,
            'Esforço cortante': False,
            'Armadura longitudinal': False,
            'Armadura transversal': False,
            'Valores': False,
        }
        grupos['Deformada'] = linhas_deformada
        # Inicialização: artistas não controlados diretamente ficam ocultos.
        for nome, artistas in grupos.items():
            visivel = estados.get(nome, False)
            for artista in artistas:
                artista.set_visible(visivel)

        margem_x = max(dx * 0.025, escala * 0.015)
        margem_y = max(dy * 0.15, escala * 0.065)
        limites = ((xmin-margem_x, xmax+margem_x),
                   (ymin-margem_y, ymax+margem_y))
        ax.set_xlim(*limites[0])
        ax.set_ylim(*limites[1])
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('x (unidade do modelo)')
        ax.set_ylabel('y (unidade do modelo)')
        ax.set_title('Modelo estrutural — geometria e apoios')
        ax.grid(False, alpha=.18)

        # Zoom com a roda do mouse, centrado na posição do cursor.
        # Funciona apenas sobre a área do gráfico, não sobre os controles.
        def zoom_scroll(event):
            if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                return
            if event.button == 'up':
                fator = 1.0 / 1.2  # Aproximar.
            elif event.button == 'down':
                fator = 1.2        # Afastar.
            else:
                return

            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            xc, yc = float(event.xdata), float(event.ydata)
            # Mantém o ponto sob o cursor na mesma posição relativa.
            novo_x0 = xc - (xc - x0) * fator
            novo_x1 = xc + (x1 - xc) * fator
            novo_y0 = yc - (yc - y0) * fator
            novo_y1 = yc + (y1 - yc) * fator
            # Evita limites degenerados em aproximações extremas.
            if abs(novo_x1 - novo_x0) < 1e-9 or abs(novo_y1 - novo_y0) < 1e-9:
                return
            ax.set_xlim(novo_x0, novo_x1)
            ax.set_ylim(novo_y0, novo_y1)
            fig.canvas.draw_idle()

        # Guarda a referência à conexão para permitir desconexão futura.
        self._zoom_scroll_cid = fig.canvas.mpl_connect('scroll_event', zoom_scroll)

        if controles:
            area = fig.add_axes([.79, .34, .20, .63])
            area.set_title('Exibir', fontsize=10, loc='left')
            nomes = [
                'Geometria da viga', 'Nós', 'Números dos nós',
                'Números das barras', 'Apoios', 'Cargas distribuídas',
                'Forças nodais', 'Momentos nodais', 'Deformada',
                'Momento fletor', 'Esforço Normal', 'Esforço cortante',
                'Armadura longitudinal', 'Armadura transversal', 'Valores',
            ]
            nome_interno = {'Esforço Normal': 'Esforço normal'}
            # Ícones quadrados de 12 pt (Matplotlib >= 3.7).
            # O parâmetro `s` corresponde à área em pontos quadrados.
            checks = CheckButtons(
                area, nomes, [estados[nome_interno.get(n, n)] for n in nomes],
                frame_props={'s': 12**2, 'linewidth': 1.4},
                check_props={'s': 12**2, 'linewidth': 1.8},
            )
            for spine in area.spines.values():
                spine.set_visible(False)

            def atualizar_valores():
                """Exibe valores somente dos diagramas ativos."""
                mostrar = estados['Valores']
                momento_ativo = mostrar and estados['Momento fletor']
                for artista in grupos['Valores dos momentos']:
                    artista.set_visible(momento_ativo)
                # Os extremos por trecho já fazem parte dos valores comuns;
                # não duplicar os rótulos globais do momento.
                for artista in grupos['Extremos do momento']:
                    artista.set_visible(False)
                for tipo, info in diagramas_nv.items():
                    ativo = mostrar and estados[tipo]
                    for artista in grupos[info['nome_valores']]:
                        artista.set_visible(ativo)
                for tipo, info in diagramas_armadura.items():
                    ativo = mostrar and estados[tipo]
                    for artista in grupos[info['nome_valores']]:
                        artista.set_visible(ativo)

            def alternar(rotulo):
                nome = nome_interno.get(rotulo, rotulo)
                estados[nome] = not estados[nome]
                if nome != 'Valores':
                    for artista in grupos[nome]:
                        artista.set_visible(estados[nome])
                # Ao ativar diagramas, oculta as cargas para evitar sobreposição.
                if (nome in ('Momento fletor', 'Esforço normal',
                             'Esforço cortante', 'Armadura longitudinal',
                             'Armadura transversal') and estados[nome]
                        and estados['Cargas distribuídas']):
                    checks.set_active(nomes.index('Cargas distribuídas'))
                atualizar_valores()
                fig.canvas.draw_idle()

            checks.on_clicked(alternar)
            area_reset = fig.add_axes([.80, .265, .17, .055])
            reset = Button(area_reset, 'Restaurar vista')

            def restaurar(_):
                ax.set_xlim(*limites[0])
                ax.set_ylim(*limites[1])
                fig.canvas.draw_idle()

            reset.on_clicked(restaurar)

            self._controles = [checks, reset]

            area_escala_m = fig.add_axes([.82, .215, .14, .020])
            slider_m = Slider(area_escala_m, 'M', 0.0, 5.0,
                              valinit=1.0, valfmt='%.2f')
            slider_m.on_changed(lambda fator: atualizar_momento(escala_momento * fator))
            self._controles.append(slider_m)
            # Sliders independentes dos esforços N e V.
            for tipo, sigla, y_slider in (
                ('Esforço normal', 'N', 0.17),
                ('Esforço cortante', 'V', 0.125),
            ):
                info = diagramas_nv[tipo]
                area_nv = fig.add_axes([.82, y_slider, .14, .020])
                controle_nv = Slider(
                    area_nv, sigla, 0.0,
                    5.0,
                    valinit=1.0, valfmt='%.2f'
                )
                controle_nv.on_changed(
                    lambda fator, atualizar=info['atualizar'], base=info['fator']:
                    atualizar(base * fator)
                )
                self._controles.append(controle_nv)



            # Uma escala para ASL/ASL2 e outra para ASW_S.
            for titulo, sigla, y_slider in (
                ('Armadura longitudinal', 'Asl', 0.100),
                ('Armadura transversal', 'Asw/s', 0.060),
            ):
                info = diagramas_armadura[titulo]
                area_arm = fig.add_axes([.82, y_slider, .14, .014])
                controle_arm = Slider(
                    area_arm, sigla, 0.0,
                    5.0,
                    valinit=1.0, valfmt='%.2f'
                )
                controle_arm.on_changed(
                    lambda fator, atualizar=info['atualizar'], base=info['fator']:
                    atualizar(base * fator)
                )
                self._controles.append(controle_arm)

            if linhas_deformada:
                # Controle de amplificação.
                area_slider = fig.add_axes([
                    0.82, 0.019, 0.14, 0.014
                ])

                slider = Slider(
                    ax=area_slider,
                    label="Def.",
                    valmin=0.0,
                    valmax=5.0,
                    valinit=1.0,
                    valfmt="%.2f"
                )

                slider.on_changed(lambda fator: atualizar_deformada(escala_auto * fator))

                self._controles.append(slider)

                # Informação sobre os deslocamentos calculados.
                no_max = max(
                    nos,
                    key=lambda no: np.hypot(
                        *deslocamento_no(no)[:2]
                    )
                )

                ux, uy, rz = deslocamento_no(no_max)

                texto_resultados = (
                    f"Nó {no_max.num}\n"
                    f"ux = {ux:.4e}\n"
                    f"uy = {uy:.4e}\n"
                    f"θz = {rz:.4e} rad"
                )

                fig.text(
                    0.08, 0.010,
                    texto_resultados,
                    fontsize=8,
                    va="bottom",
                    family="monospace"
                )


        return fig, ax


    def desenhar_deformada(self, ax=None, escala_inicial=None,
                          pontos_por_barra=30):
        """Desenha a deformada de pórtico plano com interpolação de Hermite."""
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider, CheckButtons

        if not hasattr(self.modelo, "portico"):
            raise ValueError("O modelo não possui a estrutura calculada.")

        u = np.asarray(self.modelo.portico.u, dtype=float).ravel()
        nos = self.modelo.nos
        barras = self.modelo.barras

        if len(u) != 3 * len(nos):
            raise ValueError("Dimensão incompatível do vetor de deslocamentos.")

        if not np.all(np.isfinite(u)):
            raise ValueError("O vetor de deslocamentos contém valores inválidos.")

        if ax is None:
            fig, ax = plt.subplots(figsize=(13, 5.5))
            ax.set_aspect("equal")
            ax.set_xlabel("x (unidade do modelo)")
            ax.set_ylabel("y (unidade do modelo)")
            ax.set_title("Deformada da estrutura")
            ax.grid(alpha=0.15)
            for barra in barras:
                ax.plot(
                    [barra.noi.x, barra.noj.x],
                    [barra.noi.y, barra.noj.y],
                    color="#2866a6", linewidth=1.6
                )
        else:
            fig = ax.figure

        xs = np.array([no.x for no in nos])
        ys = np.array([no.y for no in nos])
        dimensao = max(np.ptp(xs), np.ptp(ys), 1.0)
        max_u = np.max(np.hypot(u[0::3], u[1::3]))

        escala_auto = (
            0.05 * dimensao / max_u if max_u > 1e-14 else 1.0
        )
        if escala_inicial is None:
            escala_inicial = escala_auto

        escala_max = max(escala_inicial * 5.0, 1.0)
        linhas = []

        def coordenadas_deformadas(barra, fator):
            xi, yi = barra.noi.x, barra.noi.y
            xj, yj = barra.noj.x, barra.noj.y

            dx, dy = xj - xi, yj - yi
            L = np.hypot(dx, dy)

            if L <= 0:
                raise ValueError("Barra com comprimento nulo.")

            c, s = dx / L, dy / L

            def deslocamentos(no):
                i = 3 * (no.num - 1)
                return u[i:i + 3]

            ui = deslocamentos(barra.noi)
            uj = deslocamentos(barra.noj)

            # Transformação global -> local.
            ui_local = np.array([
                c * ui[0] + s * ui[1],
                -s * ui[0] + c * ui[1],
                ui[2]
            ])
            uj_local = np.array([
                c * uj[0] + s * uj[1],
                -s * uj[0] + c * uj[1],
                uj[2]
            ])

            u1, v1, t1 = ui_local
            u2, v2, t2 = uj_local

            r = np.linspace(0.0, 1.0, pontos_por_barra)

            # Interpolação axial linear.
            ua = (1.0 - r) * u1 + r * u2

            # Funções de forma de Hermite.
            H1 = 1 - 3*r**2 + 2*r**3
            H2 = L * (r - 2*r**2 + r**3)
            H3 = 3*r**2 - 2*r**3
            H4 = L * (-r**2 + r**3)

            va = H1*v1 + H2*t1 + H3*v2 + H4*t2

            # Coordenadas locais deformadas.
            xl = L*r + fator*ua
            yl = fator*va

            # Transformação local -> global.
            x = xi + c*xl - s*yl
            y = yi + s*xl + c*yl

            return x, y

        for barra in barras:
            x, y = coordenadas_deformadas(barra, escala_inicial)
            linha, = ax.plot(
                x, y,
                color="crimson",
                linewidth=2.0,
                zorder=5
            )
            linhas.append(linha)

        # Reserva espaço abaixo da figura para o controle.
        fig.subplots_adjust(bottom=0.23)

        ax_slider = fig.add_axes([0.20, 0.08, 0.57, 0.035])

        slider = Slider(
            ax=ax_slider,
            label="Amplificação",
            valmin=0.0,
            valmax=escala_max,
            valinit=escala_inicial,
            valfmt="%.2f"
        )

        def atualizar(fator):
            for barra, linha in zip(barras, linhas):
                x, y = coordenadas_deformadas(barra, fator)
                linha.set_data(x, y)

            fig.canvas.draw_idle()

        slider.on_changed(atualizar)

        # Controle para mostrar ou ocultar a deformada.
        ax_check = fig.add_axes([0.80, 0.07, 0.16, 0.07])
        check = CheckButtons(ax_check, ["Deformada"], [True])

        for spine in ax_check.spines.values():
            spine.set_visible(False)

        def alternar(_):
            visivel = not linhas[0].get_visible()
            for linha in linhas:
                linha.set_visible(visivel)
            fig.canvas.draw_idle()

        check.on_clicked(alternar)

        # Preserva referências aos controles interativos.
        if not hasattr(self, "_controles"):
            self._controles = []
        self._controles.extend([slider, check])

        return fig, ax
