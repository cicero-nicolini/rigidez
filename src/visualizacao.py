"""Visualização de modelos de pórtico plano com Matplotlib."""

# Funções legadas preservadas para compatibilidade com os cálculos.
def print_matriz(m) -> None:
    """
        Imprime uma matriz com índices e colunas alinhados.

        Args:
            m: Matriz bidimensional de valores numéricos.
    """

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
    """
        Imprime um vetor numérico com índices alinhados.

        Args:
            v: Vetor unidimensional de valores numéricos.
    """

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
        """
            Inicializa os atributos da visualização estrutural.

            Args:
                modelo: Modelo com nós, barras e resultados da análise.
        """

        self.modelo = modelo
        self._controles = []  # Mantém referências aos widgets Matplotlib.
        self.niveis_n1 = 0
        self.niveis_n2 = 0
        self.armaduras_escalonadas = []

    @staticmethod
    def _trechos_geometricos(barras):
        """
            Agrupa barras reais consecutivas, colineares e de mesma seção.

            Args:
                barras: Sequência de barras do modelo estrutural.

            Returns:
                Lista de grupos de barras geometricamente compatíveis.
        """
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
        """
            Seleciona pontos de identificação dos valores ao longo dos trechos.

            Args:
                grupos: Trechos de barras agrupadas por continuidade geométrica.
                valores_por_barra: Dicionário que associa cada barra aos valores
                    nas extremidades inicial e final.

            Returns:
                Lista de tuplas contendo nó, valor e barra de referência.
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
        """
            Obtém a interpolação gráfica dos momentos internos da barra.

            Args:
                barra: Elemento estrutural com esforços de extremidade.
                n_pontos: Número de pontos usados na representação.

            Returns:
                Tupla com as posições locais e os momentos interpolados.

            Raises:
                ValueError: Se o comprimento ou o vetor de esforços for inválido.

            Notes:
                A interpolação é linear entre -fl[2] e fl[5] e não
                representa necessariamente a distribuição exata sob cargas.
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
                       mostrar_apoios=True, ax=None, controles=True,
                       estados_iniciais=None, escalas_iniciais=None,
                       checkbox_iniciais=None, mostrar_controles_escala=True,
                       mostrar_geometria=None, mostrar_momento=None,
                       mostrar_cortante=None, mostrar_normal=None,
                       mostrar_deformada=None, mostrar_armaduras=None,
                       mostrar_valores=None, mostrar_escalas=None,
                       mostrar_modelo=None, mostrar_pontos_nos=None,
                       mostrar_numeros_nos=None, mostrar_numeros_barras=None,
                       mostrar_cargas_distribuidas=None,
                       mostrar_forcas_nodais=None, mostrar_momentos_nodais=None,
                       mostrar_diagrama_as_longitudinal=None,
                       mostrar_diagrama_as_transversal=None,
                       mostrar_armadura_longitudinal=None,
                       mostrar_armadura_transversal=None):
        """
            Desenha a estrutura e configura os controles interativos.

            Args:
                mostrar_nos: Exibe inicialmente os números dos nós.
                mostrar_barras: Exibe inicialmente os números das barras.
                mostrar_apoios: Exibe inicialmente os apoios.
                ax: Eixos Matplotlib opcionais.
                controles: Habilita os controles da interface.
                estados_iniciais: Estados de visibilidade a restaurar.
                escalas_iniciais: Escalas iniciais dos diagramas.
                checkbox_iniciais: Dicionário com os estados iniciais dos
                    checkboxes, usando os nomes apresentados na interface.
                    Os itens omitidos mantêm seus valores padrão.
                mostrar_controles_escala: Controla os sliders de amplificação.
                mostrar_geometria: Estado inicial da geometria da viga.
                mostrar_momento: Estado inicial do diagrama de momento fletor.
                mostrar_cortante: Estado inicial do diagrama de cortante.
                mostrar_normal: Estado inicial do diagrama de esforço normal.
                mostrar_deformada: Estado inicial da deformada.
                mostrar_armaduras: Estado inicial dos detalhamentos longitudinal
                    e transversal (não dos diagramas de áreas de aço).
                mostrar_valores: Estado inicial dos rótulos numéricos.
                mostrar_escalas: Exibe os sliders de escala quando True.
                    Se omitido, utiliza mostrar_controles_escala.
                mostrar_modelo: Exibe as barras do modelo estrutural.
                mostrar_pontos_nos: Exibe os pontos que representam os nós.
                mostrar_numeros_nos: Exibe a numeração dos nós. Quando omitido,
                    utiliza o valor de mostrar_nos (compatibilidade).
                mostrar_numeros_barras: Exibe a numeração das barras. Quando
                    omitido, utiliza mostrar_barras (compatibilidade).
                mostrar_cargas_distribuidas: Exibe as cargas distribuídas.
                mostrar_forcas_nodais: Exibe as forças nodais.
                mostrar_momentos_nodais: Exibe os momentos nodais.
                mostrar_diagrama_as_longitudinal: Exibe o diagrama de Asl/Asl2.
                mostrar_diagrama_as_transversal: Exibe o diagrama de Asw/s.
                mostrar_armadura_longitudinal: Exibe o detalhamento longitudinal.
                mostrar_armadura_transversal: Exibe o detalhamento transversal.

                Os parâmetros de visibilidade com valor None preservam o
                estado padrão, ou o estado restaurado da figura.

                Os parâmetros opcionais com valor None mantêm os estados
                padrão. estados_iniciais permite restaurar a interface ao
                alternar entre esforços característicos e majorados.

            Returns:
                Figura e eixos Matplotlib da visualização.

            Raises:
                ValueError: Se o modelo não contiver nós ou barras.
        """

        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.patches import Polygon, Rectangle, Circle
        from matplotlib.widgets import CheckButtons, Button, Slider, RadioButtons

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

        # ------------------------------------------------------------------
        # GEOMETRIA DA VIGA E ABERTURAS
        # ------------------------------------------------------------------
        # Silhueta única da viga e dos pilares, com uma ou várias aberturas.
        # Compatível com Modelo (abertura única) e ModeloFuros.
        geometria_viga = []
        m = self.modelo
        parametros = ('ap1', 'ap2', 'l', 'h')
        if all(hasattr(m, nome) for nome in parametros):
            from matplotlib.path import Path as MplPath
            from matplotlib.patches import PathPatch

            h = float(m.h)
            ap1, ap2 = float(m.ap1), float(m.ap2)
            x0 = ap1 / 2.0
            x1 = x0 + float(m.l)
            y0, y1 = -h / 2.0, h / 2.0

            # Os furos têm coordenadas relativas ao início do vão.
            if hasattr(m, 'furos'):
                furos = list(m.furos)
            elif all(hasattr(m, k) for k in ('lf', 'b', 'a', 'c')):
                furos = [m]
            else:
                furos = []

            if h > 0 and ap1 > 0 and ap2 > 0 and x1 > x0:
                apoio_esq = min(nos, key=lambda n: n.x)
                apoio_dir = max(nos, key=lambda n: n.x)
                xl = float(apoio_esq.x) - ap1 / 2.0
                xle = float(apoio_esq.x) + ap1 / 2.0
                xrd = float(apoio_dir.x) - ap2 / 2.0
                xr = float(apoio_dir.x) + ap2 / 2.0
                ye_inf = float(apoio_esq.y) - h
                ye_sup = float(apoio_esq.y) + h
                yd_inf = float(apoio_dir.y) - h
                yd_sup = float(apoio_dir.y) + h

                contorno = [
                    (xl, ye_inf), (xl, ye_sup),
                    (xle, ye_sup), (xle, y1),
                    (xrd, y1), (xrd, yd_sup),
                    (xr, yd_sup), (xr, yd_inf),
                    (xrd, yd_inf), (xrd, y0),
                    (xle, y0), (xle, ye_inf),
                ]

                def caminho_fechado(pontos):
                    return MplPath(
                        pontos + [pontos[0]],
                        [MplPath.MOVETO] + [MplPath.LINETO] * (len(pontos)-1)
                        + [MplPath.CLOSEPOLY]
                    )

                caminho = caminho_fechado(contorno)
                for furo in furos:
                    xc = x0 + float(furo.lf)
                    xa = xc - float(furo.b) / 2.0
                    xb = xc + float(furo.b) / 2.0
                    ya = y0 + float(furo.c)
                    yb = ya + float(furo.a)
                    if not (xle < xa < xb < xrd and y0 < ya < yb < y1):
                        continue
                    # Contorno externo é horário; os vazios são anti-horários.
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

        # ------------------------------------------------------------------
        # GRUPOS DE ARTISTAS E MODELO ESTRUTURAL
        # ------------------------------------------------------------------
        grupos = {'Modelo estrutural': [], 'Nós': [], 'Números dos nós': [], 'Números das barras': [],
                  'Apoios': [], 'Cargas distribuídas': [], 'Forças nodais': [], 'Momentos nodais': []}

        for barra in barras:
            xi, yi = barra.noi.x, barra.noi.y
            xj, yj = barra.noj.x, barra.noj.y
            linha_modelo, = ax.plot([xi, xj], [yi, yj], color='#2866a6', lw=2.1, zorder=2)
            grupos['Modelo estrutural'].append(linha_modelo)
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
            if tx and ty and rz:  # engaste: linha vertical hachurada (opção 3)
                # A hachura fica do lado externo do vão, tanto à esquerda
                # quanto à direita. Todas as linhas pertencem ao grupo Apoios.
                centro_x = (xmin + xmax) / 2.0
                lado_externo = -1 if no.x <= centro_x else 1
                meia_altura = 1.5 * tam_apoio
                linha, = ax.plot(
                    [no.x, no.x],
                    [no.y - meia_altura, no.y + meia_altura],
                    color='black', linewidth=2.2, zorder=6,
                    solid_capstyle='butt'
                )
                artistas.append(linha)
                n_hachuras = 8
                passo = 2 * meia_altura / (n_hachuras + 1)
                for i in range(1, n_hachuras + 1):
                    y0 = no.y - meia_altura + i * passo
                    # Traços diagonais apenas para fora da viga.
                    traco, = ax.plot(
                        [no.x, no.x + lado_externo * 0.55 * tam_apoio],
                        [y0, y0 - 0.55 * passo],
                        color='black', linewidth=1.1, zorder=6
                    )
                    artistas.append(traco)
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
            # Linhas do engaste já foram adicionadas por ax.plot().
            # Os demais símbolos são patches e precisam de ax.add_patch().
            from matplotlib.patches import Patch
            for artista in artistas:
                if isinstance(artista, Patch):
                    ax.add_patch(artista)
            grupos['Apoios'].extend(artistas)

        # ------------------------------------------------------------------
        # CARREGAMENTOS DISTRIBUÍDOS
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # FORÇAS E MOMENTOS NODAIS
        # ------------------------------------------------------------------
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



        # ------------------------------------------------------------------
        # DETALHAMENTO DAS ARMADURAS LONGITUDINAIS
        # ------------------------------------------------------------------
        # Detalhamento longitudinal: cores individuais e legenda fixa na área de desenho.
        grupos['Detalhamento longitudinal'] = []
        if hasattr(m, 'ap1') and hasattr(m, 'l') and hasattr(m, 'h'):
            from .armaduras_longitudinais import propor_longitudinais
            from .armaduras_banzos import propor_banzos

            propostas = propor_longitudinais(m)
            propostas_banzos = propor_banzos(m)
            artistas_banzos = []
            anotacoes_banzos = {}
            todas = [*propostas, *propostas_banzos]
            # Paleta discreta, com distinção suficiente entre as posições.
            cores = ('#c0392b', '#21618c', '#148f77', '#8e44ad',
                     '#b9770e', '#af145b', '#2874a6', '#6c7a21',
                     '#7d3c98', '#a04000', '#117864', '#566573')
            cores_por_posicao = {p.posicao: cores[i % len(cores)]
                                 for i, p in enumerate(todas)}
            cobrimento_esquematico = min(3.5, float(m.h) * 0.15)
            # Faixas de detalhamento externas à geometria: uma pista por nível.
            # Coordenadas são apenas gráficas; não alteram a posição física das barras.
            distancia_inicial = max(48.0, float(m.h) * 0.60)
            distancia_niveis = max(48.0, float(m.h) * 0.60)
            def y_faixa(face, indice):
                sinal = 1.0 if face == 'superior' else -1.0
                return sinal * (float(m.h) / 2.0 + distancia_inicial
                                + indice * distancia_niveis)

            # Os artistas de N1/N2 são atualizados sem alterar os banzos.
            artistas_n1_n2 = []

            def atualizar_faixas_n1_n2():
                from .escalonamento_n1_n2 import propor_escalonamento, continuidade_nos_furos
                for artista in artistas_n1_n2:
                    artista.remove()
                    if artista in grupos['Detalhamento longitudinal']:
                        grupos['Detalhamento longitudinal'].remove(artista)
                artistas_n1_n2.clear()
                faixas = propor_escalonamento(
                    m, niveis_n1=self.niveis_n1, niveis_n2=self.niveis_n2)
                faixas, areas_continuas = continuidade_nos_furos(m, faixas)
                self.armaduras_escalonadas = faixas
                # Redimensiona apenas os complementos dos banzos externos,
                # descontando a armadura N1/N2 que atravessa cada furo.
                novos_banzos = propor_banzos(
                    m, areas_continuas_por_furo=areas_continuas)
                for artista in artistas_banzos:
                    artista.remove()
                    if artista in grupos['Detalhamento longitudinal']:
                        grupos['Detalhamento longitudinal'].remove(artista)
                artistas_banzos.clear()
                self.armaduras_banzos = novos_banzos
                for proposta_b in novos_banzos:
                    # Armaduras dos banzos permanecem na posição original,
                    # junto à abertura, independentemente das pistas N1/N2.
                    y_b = proposta_b.y_cm
                    afast = max(3.0, min(5.0, float(m.h) * 0.07))
                    if proposta_b.banzo == 'inferior' and proposta_b.face == 'inferior':
                        y_b += afast
                    elif proposta_b.banzo == 'superior' and proposta_b.face == 'superior':
                        y_b -= afast
                    cor_b = cores_por_posicao.get(proposta_b.posicao, '#566573')
                    linha_b, = ax.plot(
                        [proposta_b.inicio_cm, proposta_b.fim_cm], [y_b, y_b],
                        color=cor_b, lw=2.6, zorder=14, visible=False)
                    rotulo_b = ax.annotate(
                        proposta_b.posicao,
                        ((proposta_b.inicio_cm + proposta_b.fim_cm) / 2, y_b),
                        xytext=(0, 10 if proposta_b.face == 'superior' else -13),
                        textcoords='offset points', ha='center', fontsize=8,
                        weight='bold', color=cor_b, visible=False, zorder=16,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=.5))
                    artistas_banzos.extend((linha_b, rotulo_b))
                    grupos['Detalhamento longitudinal'].extend((linha_b, rotulo_b))
                for pos, anot in anotacoes_banzos.items():
                    proposta_b = next((b for b in novos_banzos if b.posicao == pos), None)
                    anot.set_text(
                        f'{pos}  {proposta_b.quantidade} Ø{proposta_b.bitola_mm:g} '
                        f'C={proposta_b.comprimento_geometrico_cm:g} cm '
                        f'(As={proposta_b.as_fornecida_cm2:.2f} cm²)'
                        if proposta_b is not None else '')
                # Reserva N1 e N2 para as barras de base. As faixas
                # escalonadas recebem posições exclusivas, sem conflitar
                # com as posições já utilizadas pelas armaduras dos banzos.
                posicoes_ocupadas = {
                    proposta.posicao for proposta in novos_banzos
                } | {'N1', 'N2'}
                proximo_numero = 3

                def posicao_da_faixa(faixa):
                    """Define a identificação gráfica exclusiva de cada faixa."""
                    nonlocal proximo_numero

                    if faixa.nivel in (-1, 0):
                        return faixa.posicao

                    while f'N{proximo_numero}' in posicoes_ocupadas:
                        proximo_numero += 1

                    posicao = f'N{proximo_numero}'
                    posicoes_ocupadas.add(posicao)
                    proximo_numero += 1
                    return posicao

                for faixa in faixas:
                    posicao_grafica = posicao_da_faixa(faixa)
                    # Pista 0: construtiva. Pistas seguintes: níveis resistentes.
                    indice_pista = 0 if faixa.nivel == -1 else max(0, faixa.nivel) + 1
                    y_arm = y_faixa(faixa.face, indice_pista)
                    cor = cores_por_posicao.get(
                        faixa.posicao, '#566573'
                    )
                    linha, = ax.plot(
                        [faixa.inicio_cm, faixa.fim_cm], [y_arm, y_arm],
                        color=cor, lw=2.6, zorder=12,
                        visible=False)
                    rotulo = ax.annotate(
                        f'{faixa.quantidade} {posicao_grafica} Ø{faixa.bitola_mm:g} '
                        f'C={faixa.fim_cm - faixa.inicio_cm:.1f} cm',
                        ((faixa.inicio_cm + faixa.fim_cm) / 2, y_arm),
                        xytext=(0, 10),
                        textcoords='offset points', ha='center', va='bottom', fontsize=8,
                        weight='bold', color=cor, visible=linha.get_visible(),
                        zorder=16,
                        bbox=dict(facecolor='white', edgecolor='none', alpha=.78, pad=.5))
                    artistas_n1_n2.extend((linha, rotulo))
                    grupos['Detalhamento longitudinal'].extend((linha, rotulo))

            atualizar_faixas_n1_n2()
            self.armaduras_longitudinais = propostas

            # A relação inferior esquerda fica reservada às armaduras dos banzos.
            # As faixas externas já exibem posição, quantidade, bitola e comprimento.
            for i, proposta in enumerate(propostas_banzos):
                texto = (f'{proposta.posicao}  {proposta.quantidade} '
                         f'Ø{proposta.bitola_mm:g}  '
                         f'C={proposta.comprimento_geometrico_cm:g} cm '
                         f'(As={proposta.as_fornecida_cm2:.2f} cm²)')
                anot = ax.annotate(
                    texto, xy=(0.018, 0.035), xycoords='axes fraction',
                    xytext=(0, 12 * (len(propostas_banzos) - 1 - i)),
                    textcoords='offset points', ha='left', va='bottom',
                    color=cores_por_posicao[proposta.posicao],
                    fontsize=8, family='monospace', weight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=1.0),
                    visible=False, zorder=30, annotation_clip=False)
                grupos['Detalhamento longitudinal'].append(anot)
                anotacoes_banzos[proposta.posicao] = anot

        # ------------------------------------------------------------------
        # DETALHAMENTO DAS ARMADURAS TRANSVERSAIS
        # ------------------------------------------------------------------
        # Estribos em seção cheia e nos banzos, com relação à direita.
        grupos['Detalhamento transversal'] = []
        if hasattr(m, 'ap1') and hasattr(m, 'l') and (hasattr(m, 'furos') or hasattr(m, 'lf')):
            from .armaduras_transversais import propor_transversais
            from .armaduras_estribos_banzos import propor_estribos_banzos
            propostas_estribos = propor_transversais(m)
            propostas_banzos_estribos = propor_estribos_banzos(
                m, primeira_posicao=len(propostas_estribos) + 1)
            todos_estribos = [*propostas_estribos, *propostas_banzos_estribos]
            h = float(m.h)
            cobrimento = min(3.5, h * 0.15)
            topo, base = h / 2 - cobrimento, -h / 2 + cobrimento
            cores_estribos = ('#8050a0', '#b05b20', '#167f76', '#b03267',
                              '#4169a2', '#737c18', '#a04a36', '#6b55a8')
            for i, proposta in enumerate(todos_estribos):
                cor = cores_estribos[i % len(cores_estribos)]
                if hasattr(proposta, 'y_inferior_cm'):
                    ybaixo, yalto = proposta.y_inferior_cm, proposta.y_superior_cm
                else:
                    ybaixo, yalto = base, topo
                xs = np.linspace(proposta.inicio_cm, proposta.fim_cm,
                                 proposta.quantidade)
                for x in xs:
                    linha, = ax.plot([x, x], [ybaixo, yalto], color=cor,
                                     lw=.85, zorder=11, visible=False)
                    grupos['Detalhamento transversal'].append(linha)
                centro = (proposta.inicio_cm + proposta.fim_cm) / 2
                rotulo = ax.annotate(
                    proposta.posicao, (centro, (ybaixo + yalto) / 2),
                    ha='center', va='center', fontsize=8, weight='bold',
                    color=cor, zorder=17, visible=False,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=1))
                grupos['Detalhamento transversal'].append(rotulo)
                anot = ax.annotate(
                    (f'{proposta.posicao} Ø{proposta.bitola_mm:g} '
                     f'c/{proposta.espacamento_cm:g} '
                     f'(Asws={proposta.asw_fornecida_cm2_m:.2f} cm²/m)'),
                    xy=(.982, .035), xycoords='axes fraction',
                    xytext=(0, 12 * (len(todos_estribos) - 1 - i)),
                    textcoords='offset points', ha='right', va='bottom',
                    color=cor, fontsize=8, family='monospace', weight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=1),
                    visible=False, zorder=30, annotation_clip=False)
                grupos['Detalhamento transversal'].append(anot)
            self.armaduras_transversais = propostas_estribos
            self.armaduras_estribos_banzos = propostas_banzos_estribos

        grupos['Geometria da viga'] = geometria_viga

        # ------------------------------------------------------------------
        # CONFIGURAÇÃO DE VISIBILIDADE E CONTROLES
        # ------------------------------------------------------------------
        # Os rótulos numéricos são controlados por uma única opção "Valores".
        # As coleções de artistas continuam separadas internamente.
        estados = {
            'Geometria da viga': False,
            'Modelo estrutural': True,
            'Nós': True,
            'Números dos nós': mostrar_nos,
            'Números das barras': mostrar_barras,
            'Apoios': mostrar_apoios if estados_iniciais is None else None,
            'Cargas distribuídas': True,
            'Forças nodais': True,
            'Momentos nodais': True,
            'Deformada': False,
            'Momento fletor': False,
            'Esforço normal': False,
            'Esforço cortante': False,
            'Armadura longitudinal': False,
            'Detalhamento longitudinal': False,
            'Detalhamento transversal': False,
            'Armadura transversal': False,
            'Valores': False,
        }
        # Relaciona os rótulos da interface aos grupos de artistas internos.
        # Os dois diagramas de armadura têm nomes distintos do detalhamento.
        nomes_internos_checkbox = {
            'Esforço Normal': 'Esforço normal',
            'Diagrama As longitudinal': 'Armadura longitudinal',
            'Diagrama As transversal': 'Armadura transversal',
            'Armadura longitudinal': 'Detalhamento longitudinal',
            'Armadura transversal': 'Detalhamento transversal',
        }

        # Ao reconstruir a figura, preserva os estados sem acionar callbacks.
        if estados_iniciais is not None:
            estados.update({chave: bool(valor) for chave, valor in
                            estados_iniciais.items() if chave in estados})

        # Os estados fornecidos explicitamente têm prioridade sobre os padrões.
        # A configuração antiga continua disponível para compatibilidade.
        if checkbox_iniciais is not None:
            if not isinstance(checkbox_iniciais, dict):
                raise TypeError('checkbox_iniciais deve ser um dicionário.')

            for rotulo, marcado in checkbox_iniciais.items():
                nome = nomes_internos_checkbox.get(rotulo, rotulo)
                if nome not in estados:
                    raise ValueError(f'Checkbox desconhecido: {rotulo!r}.')
                estados[nome] = bool(marcado)

        # Configuração simplificada: parâmetros booleanos da chamada.
        # None significa que o estado padrão (ou restaurado) é preservado.
        opcoes_simples = {
            # Geometria, modelo e identificação dos elementos.
            'Geometria da viga': mostrar_geometria,
            'Modelo estrutural': mostrar_modelo,
            'Nós': mostrar_pontos_nos,
            'Números dos nós': mostrar_numeros_nos,
            'Números das barras': mostrar_numeros_barras,
            'Apoios': mostrar_apoios,

            # Carregamentos aplicados ao modelo.
            'Cargas distribuídas': mostrar_cargas_distribuidas,
            'Forças nodais': mostrar_forcas_nodais,
            'Momentos nodais': mostrar_momentos_nodais,

            # Resultados da análise estrutural.
            'Deformada': mostrar_deformada,
            'Momento fletor': mostrar_momento,
            'Esforço normal': mostrar_normal,
            'Esforço cortante': mostrar_cortante,

            # Diagramas de armadura e detalhamentos.
            'Armadura longitudinal': mostrar_diagrama_as_longitudinal,
            'Armadura transversal': mostrar_diagrama_as_transversal,
            'Detalhamento longitudinal': mostrar_armadura_longitudinal,
            'Detalhamento transversal': mostrar_armadura_transversal,
            'Valores': mostrar_valores,
        }

        for nome, marcado in opcoes_simples.items():
            if marcado is not None:
                estados[nome] = bool(marcado)

        if mostrar_armaduras is not None:
            estados['Detalhamento longitudinal'] = bool(mostrar_armaduras)
            estados['Detalhamento transversal'] = bool(mostrar_armaduras)

        # O nome curto tem prioridade quando os dois forem informados.
        if mostrar_escalas is not None:
            mostrar_controles_escala = bool(mostrar_escalas)

        grupos['Deformada'] = linhas_deformada
        # Inicialização: artistas não controlados diretamente ficam ocultos.
        for nome, artistas in grupos.items():
            visivel = estados.get(nome, False)
            for artista in artistas:
                artista.set_visible(visivel)

        margem_x = max(dx * 0.025, escala * 0.015)
        margem_y = max(dy * 0.15, escala * 0.065)
        # Reserva espaço para as pistas externas, inclusive após alternar níveis.
        if hasattr(m, 'h') and hasattr(m, 'ap1'):
            margem_pistas = max(38.0, float(m.h) * 0.48) + 5 * max(32.0, float(m.h) * 0.40)
            ymin = min(ymin, -float(m.h) / 2 - margem_pistas)
            ymax = max(ymax, float(m.h) / 2 + margem_pistas)
        limites = ((xmin-margem_x, xmax+margem_x),
                   (ymin-margem_y, ymax+margem_y))
        ax.set_xlim(*limites[0])
        ax.set_ylim(*limites[1])
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('x (unidade do modelo)')
        ax.set_ylabel('y (unidade do modelo)')
        ax.set_title('Modelo estrutural — geometria e apoios')
        ax.grid(False)

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
                'Geometria da viga', 'Modelo estrutural', 'Nós', 'Números dos nós',
                'Números das barras', 'Apoios', 'Cargas distribuídas',
                'Forças nodais', 'Momentos nodais', 'Deformada',
                'Momento fletor', 'Esforço Normal', 'Esforço cortante',
                'Diagrama As longitudinal', 'Diagrama As transversal',
                'Armadura longitudinal', 'Armadura transversal', 'Valores',
            ]
            nome_interno = nomes_internos_checkbox
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
            # Os artistas de rótulos são criados ocultos. Ao reconstruir a
            # figura após trocar a combinação de esforços, restaura também
            # sua visibilidade conforme os checkboxes preservados, sem
            # precisar clicar novamente em "Valores".
            atualizar_valores()
            area_reset = fig.add_axes([.80, .265, .17, .055])
            reset = Button(area_reset, 'Restaurar vista')

            def restaurar(_):
                ax.set_xlim(*limites[0])
                ax.set_ylim(*limites[1])
                fig.canvas.draw_idle()

            reset.on_clicked(restaurar)

            self._controles = [checks, reset]

            # Controle independente de níveis em uma pequena janela.
            area_escal = fig.add_axes([.80, .306, .17, .030])
            botao_escal = Button(area_escal, 'Escalonar N1 / N2')

            def abrir_escalonamento(_):
                from matplotlib.widgets import Slider as SliderEsc
                janela = plt.figure(figsize=(7.4, 3.0))
                janela.suptitle('Escalonamento teórico — N1 e N2 (sem ancoragens)')
                a1 = janela.add_axes([.23, .59, .62, .11])
                a2 = janela.add_axes([.23, .32, .62, .11])
                controle_n1 = SliderEsc(a1, 'N1 inferior', 0, 3,
                                        valinit=self.niveis_n1, valstep=1)
                controle_n2 = SliderEsc(a2, 'N2 superior', 0, 3,
                                        valinit=self.niveis_n2, valstep=1)
                aviso = janela.text(.08, .09,
                    '0 = contínua; 1–3 = máximo de grupos adicionais.\n'
                    'Comprimentos teóricos, não utilizar para corte de barras.',
                    fontsize=9)

                def alterar(_):
                    self.niveis_n1 = int(controle_n1.val)
                    self.niveis_n2 = int(controle_n2.val)
                    atualizar_faixas_n1_n2()
                    visivel = estados['Detalhamento longitudinal']
                    for artista in (*artistas_n1_n2, *artistas_banzos):
                        artista.set_visible(visivel)
                    fig.canvas.draw_idle()

                controle_n1.on_changed(alterar)
                controle_n2.on_changed(alterar)
                # Referências mantidas para evitar coleta dos widgets.
                janela._controles_escalonamento = (controle_n1, controle_n2, aviso)
                janela.show()

            botao_escal.on_clicked(abrir_escalonamento)
            self._controles.append(botao_escal)

            # Seletor disponível somente para modelos calculados nas duas
            # combinações (opção 2). Os detalhamentos permanecem de ELU.
            if (hasattr(self.modelo, 'resultados_caracteristicos')
                    and hasattr(self.modelo, 'resultados_calculo')
                    and callable(getattr(self.modelo, 'selecionar_esforcos', None))):
                # Posiciona o seletor no canto superior direito da área
                # de plotagem, sem sobrepor os controles de escala laterais.
                area_tipo = ax.inset_axes([.785, .815, .205, .145],
                                         transform=ax.transAxes, zorder=20)
                area_tipo.set_facecolor((1, 1, 1, 0.92))
                area_tipo.set_title('Esforços exibidos', fontsize=9, loc='left')
                rotulos_tipo = ['Característicos', 'Majorados']
                tipo_atual = getattr(self.modelo, 'tipo_esforcos', 'calculo')
                seletor = RadioButtons(area_tipo, rotulos_tipo,
                                       active=0 if tipo_atual == 'caracteristicos' else 1)
                for spine in area_tipo.spines.values():
                    spine.set_visible(False)

                def trocar_tipo(rotulo):
                    novo_tipo = ('caracteristicos' if rotulo == 'Característicos'
                                 else 'calculo')
                    if novo_tipo == getattr(self.modelo, 'tipo_esforcos', 'calculo'):
                        return
                    # Reconstrói os diagramas: os artistas anteriores contêm
                    # valores numéricos calculados na criação da figura.
                    # Captura os estados internos, incluindo nomes de grupos
                    # diferentes dos rótulos apresentados nos checkboxes.
                    estados_anteriores = estados.copy()
                    escalas_anteriores = {
                        controle.label.get_text(): controle.val
                        for controle in self._controles
                        if isinstance(controle, Slider)
                    }
                    limites_atuais = (ax.get_xlim(), ax.get_ylim())
                    self.modelo.selecionar_esforcos(novo_tipo)
                    if hasattr(self, '_zoom_scroll_cid'):
                        fig.canvas.mpl_disconnect(self._zoom_scroll_cid)
                    for controle_antigo in self._controles:
                        if hasattr(controle_antigo, 'disconnect_events'):
                            controle_antigo.disconnect_events()
                    fig.clear()
                    novo_ax = fig.add_subplot(111)
                    self.desenhar_geometria(
                        ax=novo_ax, controles=controles,
                        estados_iniciais=estados_anteriores,
                        escalas_iniciais=escalas_anteriores,
                        mostrar_controles_escala=mostrar_controles_escala,
                    )
                    novo_ax.set_xlim(*limites_atuais[0])
                    novo_ax.set_ylim(*limites_atuais[1])
                    fig.canvas.draw_idle()

                seletor.on_clicked(trocar_tipo)
                self._controles.append(seletor)

            if mostrar_controles_escala:
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

            # Reaplica as escalas selecionadas pelo usuário após recriar
            # os artistas. Não altera os estados dos checkboxes.
            if escalas_iniciais:
                for controle in self._controles:
                    if isinstance(controle, Slider):
                        rotulo = controle.label.get_text()
                        if rotulo in escalas_iniciais:
                            valor = float(escalas_iniciais[rotulo])
                            if controle.valmin <= valor <= controle.valmax:
                                controle.set_val(valor)

        return fig, ax


    def desenhar_deformada(self, ax=None, escala_inicial=None,
                          pontos_por_barra=30):
        """
            Apresenta a deformada da estrutura em uma figura Matplotlib.

            Args:
                ax: Eixos Matplotlib opcionais.
                escala_inicial: Escala inicial da deformada.
                pontos_por_barra: Número de pontos de interpolação por barra.
        """
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
