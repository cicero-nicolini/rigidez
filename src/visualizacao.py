"""Visualização de modelos de pórtico plano com Matplotlib."""

SIMBOLO_BITOLA = "\u00D8"  # Ø: símbolo único para bitolas em toda a interface.

# Funções legadas preservadas para compatibilidade com os cálculos.
from .utilidades import print_matriz, print_vetor
from .detalhamento import Detalhamento
from .escalonamento_n1_n2 import calcula_momentos_decalados


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
        self.decalagem_cm = 0.0  # Decalagem geométrica da demanda (cm).
        self._decalagem_manual = False

    def desenhar_analise_estrutural(self, mostrar=True, carregamentos=()):
        """Desenha pórtico, cargas e diagramas em um único eixo interativo.

        Args:
            mostrar: Abre a janela Matplotlib quando verdadeiro.
            carregamentos: Cargas aplicadas às barras da estrutura.

        Returns:
            Tupla com figura e eixo único de visualização.

        Raises:
            ValueError: Quando faltam resultados de análise válidos.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Slider, Button
        from matplotlib.patches import FancyArrowPatch, Polygon, Circle, Arc

        estrutura = self.modelo
        nos, barras = estrutura.nos, estrutura.barras
        u = np.asarray(estrutura.u, dtype=float).ravel()
        if u.size != 3 * len(nos) or not np.all(np.isfinite(u)):
            raise ValueError("Deslocamentos inválidos ou não calculados.")
        if any(not hasattr(b, "fl") for b in barras):
            raise ValueError("Calcule os esforços das barras antes de desenhar.")

        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_axes([0.055, 0.09, 0.70, 0.83])
        ax.set_aspect("equal", adjustable="datalim")
        # Sem título na área de plotagem do Exemplo 1.
        ax.set_xlabel("x (unidade do modelo)")
        ax.set_ylabel("y (unidade do modelo)")
        ax.grid(alpha=0.13)
        extensao = max(np.ptp([n.x for n in nos]), np.ptp([n.y for n in nos]), 1.0)
        grupos = {k: [] for k in ("Geometria do modelo", "Vinculações", "Carregamentos", "Reações",
                                  "Deformada", "Esforço Normal", "Esforço Cortante", "Esforço Fletor")}
        escalas = {k: 1.0 for k in ("Carregamentos", "Reações", "Deformada", "Esforço Normal", "Esforço Cortante", "Esforço Fletor")}
        elementos_escalaveis = {k: [] for k in escalas}
        setas_escalaveis = {"Carregamentos": [], "Reações": []}
        cor_carga = "#b83280"
        cor_reacao = "#1846d0"

        def registrar(grupo, artista):
            grupos[grupo].append(artista)
            return artista

        def linha_escalavel(grupo, base, vetor, cor, largura=1.8, texto=None):
            """Registra linha e eventual rótulo para atualizar sua escala."""
            pts = base + escalas[grupo] * vetor
            art, = ax.plot(pts[:, 0], pts[:, 1], color=cor, lw=largura, zorder=4)
            registrar(grupo, art)
            rotulo = None
            if texto is not None:
                rotulo = registrar(grupo, ax.annotate(texto, pts[-1], fontsize=8,
                                                      xytext=(4, 4), textcoords="offset points",
                                                      color=cor))
            elementos_escalaveis[grupo].append((art, rotulo, base, vetor))

        def seta_carga(pos, direcao, comprimento, texto=None, grupo="Carregamentos", cor=None):
            """Desenha flecha com ponta no ponto de aplicação, escalável."""
            pos = np.asarray(pos, float)
            direcao = np.asarray(direcao, float)
            norma = np.linalg.norm(direcao)
            if norma < 1e-12:
                return
            cor = cor or (cor_reacao if grupo == "Reações" else cor_carga)
            vetor = comprimento * direcao / norma
            inicio = pos if grupo == "Reações" else pos - escalas[grupo] * vetor
            fim = pos + escalas[grupo] * vetor if grupo == "Reações" else pos
            seta = FancyArrowPatch(inicio, fim, arrowstyle='-|>',
                                    mutation_scale=12, lw=1.4,
                                    color=cor, zorder=8)
            ax.add_patch(seta)
            registrar(grupo, seta)
            rotulo = None
            if texto:
                rotulo = registrar(grupo, ax.annotate(
                    texto, fim if grupo == "Reações" else inicio, xytext=(3, 3), textcoords='offset points',
                    fontsize=8, color=cor))
            setas_escalaveis[grupo].append((seta, rotulo, pos, vetor))

        momentos_escalaveis = {"Carregamentos": [], "Reações": []}

        def desenhar_momento_escalavel(pos, valor, grupo, cor, texto):
            """Arco e ponta de seta solidários à escala do grupo."""
            pos = np.asarray(pos, dtype=float)
            raio = (.055 if grupo == "Carregamentos" else .042) * extensao
            sentido = 1 if valor > 0 else -1
            angulos = np.deg2rad(np.linspace(35, 325, 65))
            if sentido < 0:
                angulos = angulos[::-1]
            curva = np.column_stack((np.cos(angulos), np.sin(angulos))) * raio
            pontos = pos + curva
            linha, = ax.plot(pontos[:, 0], pontos[:, 1], color=cor, lw=1.7, zorder=9)
            registrar(grupo, linha)
            tangente = sentido * np.array([-np.sin(angulos[-1]), np.cos(angulos[-1])])
            vetor_ponta = curva[-1]
            ponta = pos + vetor_ponta
            seta = FancyArrowPatch(ponta - .35 * raio * tangente, ponta,
                                   arrowstyle="-|>", mutation_scale=15,
                                   lw=1.7, color=cor, zorder=10)
            ax.add_patch(seta)
            registrar(grupo, seta)
            # O rótulo acompanha a extremidade superior do arco ao escalar.
            indice_rotulo = int(np.argmax(curva[:, 1]))
            ancora_rotulo = pos + escalas[grupo] * curva[indice_rotulo]
            rotulo = registrar(grupo, ax.annotate(
                texto, ancora_rotulo, xytext=(4, 5), textcoords="offset points",
                fontsize=8, color=cor))
            momentos_escalaveis[grupo].append((linha, seta, pos, curva,
                                                tangente, raio, rotulo, indice_rotulo))

        # Geometria, numeração e diagramas de esforços.
        cores = {"Esforço Normal": "#7c3aed", "Esforço Cortante": "#e07b18",
                 "Esforço Fletor": "#16803c"}
        for barra, x, N, V, M in estrutura.calcula_solicitacoes_internas(carregamentos, 201):
            p = np.array([barra.noi.x, barra.noi.y], float)
            q = np.array([barra.noj.x, barra.noj.y], float)
            d = q-p
            L = float(np.linalg.norm(d))
            if L <= 0:
                raise ValueError(f"Barra {barra.num} possui comprimento nulo.")
            t, n = d/L, np.array([-d[1], d[0]])/L
            registrar("Geometria do modelo", ax.plot([p[0], q[0]], [p[1], q[1]],
                                                    color="0.35", lw=1.8, zorder=3)[0])
            meio = (p+q)/2
            registrar("Geometria do modelo", ax.annotate(f"B{int(barra.num)}", meio,
                        xytext=(5, 6), textcoords="offset points", fontsize=9, color="0.25"))
            dados = {"N": N, "V": V, "M": M}
            eixo = p + x[:, None] * t
            for grupo, campo in (("Esforço Normal", "N"), ("Esforço Cortante", "V"),
                                 ("Esforço Fletor", "M")):
                valores = np.asarray(dados[campo], float)
                if campo == "V":
                    valores = -valores  # Convenção de exibição do FTOOL.
                lado = -1 if campo == "M" else 1
                amp = max(np.max(np.abs(valores)), 1e-12)
                vetor = lado * 0.14 * extensao * (valores/amp)[:, None] * n
                linha_escalavel(grupo, eixo, vetor, cores[grupo], 1.8)
                # Preenchimento translúcido, conforme os diagramas da
                # visualização das vigas com aberturas.
                contorno = np.vstack((eixo, (eixo + vetor)[::-1]))
                preenchimento = Polygon(contorno, closed=True,
                                        facecolor=cores[grupo], edgecolor='none',
                                        alpha=0.13, zorder=2)
                ax.add_patch(preenchimento)
                registrar(grupo, preenchimento)
                elementos_escalaveis[grupo].append(
                    (preenchimento, None, eixo, vetor)
                )
                indices = {0, len(valores)-1}
                if campo == "M":
                    indices.update((int(np.argmin(valores)), int(np.argmax(valores))))
                    for i in range(1, len(valores)-1):
                        if (valores[i]-valores[i-1])*(valores[i+1]-valores[i]) < 0:
                            indices.add(i)
                for i in sorted(indices):
                    linha_escalavel(grupo, np.stack((eixo[i], eixo[i])),
                                   np.stack((np.zeros(2), vetor[i])), cores[grupo], .7,
                                   f"{valores[i]:.3g}")

        for no in nos:
            xy = (no.x, no.y)
            registrar("Geometria do modelo", ax.plot(*xy, "ko", ms=4, zorder=7)[0])
            registrar("Geometria do modelo", ax.annotate(f"{no.num}", xy,
                       xytext=(5, 7), textcoords="offset points", fontsize=9))
            # Vinculações com símbolos tradicionais, como no modelo de vigas.
            tx, ty, rz = bool(no.Tx), bool(no.Ty), bool(no.Rz)
            tamanho = 0.025 * extensao
            if tx and ty and rz:
                centro = (min(n.x for n in nos) + max(n.x for n in nos)) / 2
                lado = -1 if no.x <= centro else 1
                registrar("Vinculações", ax.plot(
                    [no.x, no.x], [no.y-1.4*tamanho, no.y+1.4*tamanho],
                    color='black', lw=2.2, zorder=9)[0])
                for k in range(8):
                    yy = no.y + (k-3.5)*0.34*tamanho
                    registrar("Vinculações", ax.plot(
                        [no.x, no.x+lado*.55*tamanho],
                        [yy, yy-.22*tamanho], color='black', lw=1, zorder=9)[0])
            elif tx and ty:
                registrar("Vinculações", ax.add_patch(Polygon(
                    [(no.x,no.y),(no.x-tamanho,no.y-tamanho),
                     (no.x+tamanho,no.y-tamanho)], closed=True,
                    facecolor='#dddddd', edgecolor='black', zorder=9)))
            elif ty:
                registrar("Vinculações", ax.add_patch(Polygon(
                    [(no.x,no.y),(no.x-tamanho,no.y-.8*tamanho),
                     (no.x+tamanho,no.y-.8*tamanho)], closed=True,
                    facecolor='#dddddd', edgecolor='black', zorder=9)))
                for lado in (-.5,.5):
                    registrar("Vinculações", ax.add_patch(Circle(
                        (no.x+lado*tamanho,no.y-tamanho), .16*tamanho,
                        facecolor='white', edgecolor='black', zorder=9)))
            elif tx:
                registrar("Vinculações", ax.add_patch(Polygon(
                    [(no.x,no.y),(no.x-tamanho,no.y-tamanho),
                     (no.x-tamanho,no.y+tamanho)], closed=True,
                    facecolor='#dddddd', edgecolor='black', zorder=9)))

            # Cargas nodais em coordenadas globais.
            for componente, direcao in (("Fx", np.array([1., 0.])),
                                        ("Fy", np.array([0., 1.]))):
                valor = float(getattr(no, componente, 0.0))
                if abs(valor) > 1e-12:
                    seta_carga(xy, np.sign(valor)*direcao,
                               .13*extensao, f"{componente}={valor:g}")
            momento = float(getattr(no, "Mz", 0.0))
            if abs(momento) > 1e-12:
                desenhar_momento_escalavel(xy, momento, "Carregamentos",
                                           cor_carga, f"Mz={momento:g}")

        # Reações nos graus de liberdade vinculados (vetor global Estrutura.R).
        # O método calcula_reacoes() deve ter sido executado antes da visualização.
        reacoes = np.asarray(getattr(estrutura, "R", []), dtype=float).ravel()
        if reacoes.size == 3 * len(nos):
            for no in nos:
                indice = 3 * (int(no.num) - 1)
                if indice < 0 or indice + 2 >= reacoes.size:
                    continue
                pos = np.array([no.x, no.y], dtype=float)
                for restrito, valor, direcao, nome in (
                    (no.Tx, reacoes[indice], np.array([1., 0.]), "Rx"),
                    (no.Ty, reacoes[indice + 1], np.array([0., 1.]), "Ry"),
                ):
                    if restrito and abs(valor) > 1e-10:
                        seta_carga(pos, np.sign(valor) * direcao,
                                   .11 * extensao, f"{nome}={valor:.3g}",
                                   grupo="Reações", cor=cor_reacao)
                momento_r = reacoes[indice + 2]
                if no.Rz and abs(momento_r) > 1e-10:
                    desenhar_momento_escalavel(pos, momento_r, "Reações",
                                               cor_reacao, f"Rz={momento_r:.3g}")

        # Cargas distribuídas: flechas com pontas na barra e linha de origem.
        for carga in carregamentos:
            barra = carga.barra
            p = np.array([barra.noi.x, barra.noi.y], float)
            q = np.array([barra.noj.x, barra.noj.y], float)
            L = np.linalg.norm(q-p)
            t = (q-p)/L
            n = np.array([-t[1], t[0]])
            if hasattr(carga, "lw"):
                xs = np.linspace(float(carga.a), float(carga.a+carga.lw),
                                 max(3, int(carga.lw / max(extensao*.035, 1e-9))+1))
                valores = [float(carga.w1)+(float(carga.w2)-float(carga.w1))*
                           (x-carga.a)/carga.lw for x in xs]
                maior = max(max(abs(v) for v in valores), 1e-12)
                origens = []
                for x, w in zip(xs, valores):
                    pos = p+x*t
                    direcao = -np.sign(w)*n
                    comprimento = .055*extensao*abs(w)/maior
                    if comprimento > 1e-10:
                        seta_carga(pos, direcao, comprimento)
                        origens.append((pos, direcao, comprimento))
                # Linha superior que liga as origens das flechas.
                if len(origens) >= 2:
                    posicoes = np.array([pos-compr*direcao for pos,direcao,compr in origens])
                    base = np.array([pos for pos,_,_ in origens])
                    linha_escalavel("Carregamentos", base, posicoes-base,
                                    cor_carga, .9)
                meio = p+(float(carga.a)+float(carga.lw)/2)*t
                # Posicionar o valor na linha de origem das flechas, não
                # em uma distância fixa em pontos da geometria.
                if origens:
                    indice_meio = len(origens) // 2
                    pos_m, direcao_m, comprimento_m = origens[indice_meio]
                    base_rotulo = np.asarray(pos_m, float)
                    vetor_rotulo = -comprimento_m * np.asarray(direcao_m, float)
                    pos_rotulo = base_rotulo + escalas["Carregamentos"] * vetor_rotulo
                    texto_rotulo = (f"{float(carga.w1):g}"
                                    if abs(carga.w1-carga.w2)<1e-12
                                    else f"{carga.w1:g}–{carga.w2:g}")
                    rotulo = registrar("Carregamentos", ax.annotate(
                        texto_rotulo, pos_rotulo,
                        xytext=(0, 5), textcoords='offset points',
                        color=cor_carga, fontsize=8, ha='center'))
                    elementos_escalaveis["Carregamentos"].append(
                        (None, rotulo, base_rotulo[None, :], vetor_rotulo[None, :]))
            else:
                pos = p+float(carga.a)*t
                # Componentes locais: Py positivo para baixo no exemplo.
                for valor, direcao, nome in ((float(carga.Px), t, 'Px'),
                                             (float(carga.Py), -n, 'Py')):
                    if abs(valor) > 1e-12:
                        seta_carga(pos, np.sign(valor)*direcao,
                                   .14*extensao, f"{nome}={valor:g}")

        # Deformada obtida da análise estrutural.
        max_deslocamento = float(np.max(np.hypot(u[0::3], u[1::3])))
        fator_deformada = (.06 * extensao / max_deslocamento
                           if max_deslocamento > 1e-14 else 1.0)
        deformadas = []
        for barra, original, desloc in estrutura.calcula_deformada(81):
            pontos = original + fator_deformada*desloc
            linha, = ax.plot(pontos[:, 0], pontos[:, 1], color="#1769aa",
                             lw=2, zorder=6)
            registrar("Deformada", linha)
            deformadas.append((linha, original, desloc))

        # Destaque do nó selecionado, independente dos símbolos de vinculação.
        marcador_no_selecionado, = ax.plot([], [], "o", ms=10,
            markerfacecolor="red", markeredgecolor="darkred",
            markeredgewidth=1.5, zorder=30, visible=False)

        # Clique próximo a um nó: deslocamentos no canto inferior esquerdo.
        painel_no = ax.text(.015, .025, "Clique em um nó para ver seus deslocamentos",
                            transform=ax.transAxes, fontsize=9, va="bottom",
                            ha="left", color="0.2",
                            bbox=dict(facecolor="white", alpha=.85, edgecolor="0.75"),
                            zorder=20)

        # Painel único de opções à direita.
        nomes = list(grupos)
        check_ax = fig.add_axes([0.79, 0.57, 0.20, 0.34])
        checks = CheckButtons(
            check_ax, nomes + ["Grid"],
            [True, True, True, False, False, False, False, False, True],
            frame_props={"s": 12**2, "linewidth": 1.4},
            check_props={"s": 12**2, "linewidth": 1.8},
        )
        check_ax.set_title("Exibir", fontsize=10, loc="left")
        for spine in check_ax.spines.values():
            spine.set_visible(False)
        for nome in nomes[3:]:
            for artista in grupos[nome]:
                artista.set_visible(False)

        def alternar(nome):
            if nome == "Grid":
                ax.grid(checks.get_status()[-1])
                fig.canvas.draw_idle()
                return
            itens = grupos[nome]
            if itens:
                novo = not itens[0].get_visible()
                for item in itens:
                    item.set_visible(novo)
            fig.canvas.draw_idle()
        checks.on_clicked(alternar)

        sliders = []
        rotulos_escala = {"Carregamentos": "carregamentos",
                         "Reações": "reações",
                         "Deformada": "deformada",
                         "Esforço Normal": "esforço normal",
                         "Esforço Cortante": "esforço cortante",
                         "Esforço Fletor": "momento fletor"}
        for j, grupo in enumerate(escalas):
            y_slider = 0.475 - j * 0.060
            fig.text(0.795, y_slider + 0.037, rotulos_escala[grupo],
                     fontsize=9, ha="left", va="center")
            sx = fig.add_axes([0.805, y_slider, 0.15, 0.020])
            controle = Slider(sx, "", 0.0, 5.0, valinit=1.0, valfmt="%.2f")
            def atualizar(valor, grupo=grupo):
                escalas[grupo] = valor
                for art, rotulo, base, vetor in elementos_escalaveis[grupo]:
                    pts = base + valor * vetor
                    if isinstance(art, Polygon):
                        art.set_xy(np.vstack((base, pts[::-1])))
                    elif art is not None:
                        art.set_data(pts[:, 0], pts[:, 1])
                    if rotulo is not None:
                        rotulo.xy = tuple(pts[-1])
                if grupo in setas_escalaveis:
                    for seta, rotulo, pos, vetor in setas_escalaveis[grupo]:
                        inicio = pos if grupo == "Reações" else pos - valor * vetor
                        fim = pos + valor * vetor if grupo == "Reações" else pos
                        seta.set_positions(inicio, fim)
                        if rotulo is not None:
                            rotulo.xy = tuple(fim if grupo == "Reações" else inicio)
                if grupo in momentos_escalaveis:
                    for linha, seta, centro, curva, tangente, raio, rotulo, indice_rotulo in momentos_escalaveis[grupo]:
                        pontos = centro + valor * curva
                        linha.set_data(pontos[:, 0], pontos[:, 1])
                        ponta = pontos[-1]
                        seta.set_positions(ponta - valor*.35*raio*tangente, ponta)
                        rotulo.xy = tuple(pontos[indice_rotulo])
                if grupo == "Deformada":
                    for linha, original, desloc in deformadas:
                        pontos = original + valor*fator_deformada*desloc
                        linha.set_data(pontos[:, 0], pontos[:, 1])
                fig.canvas.draw_idle()
            controle.on_changed(atualizar)
            sliders.append(controle)

        # Escala independente dos textos da plotagem, inclusive dos ocultos.
        # Preserva as diferenças entre os tamanhos originais das fontes.
        textos_plotagem = [(texto, texto.get_fontsize()) for texto in ax.texts]
        textos_plotagem.extend((texto, texto.get_fontsize()) for texto in
                               (ax.xaxis.label, ax.yaxis.label))
        tamanho_ticks_x = ax.xaxis.get_major_ticks()[0].label1.get_fontsize()
        tamanho_ticks_y = ax.yaxis.get_major_ticks()[0].label1.get_fontsize()
        fig.text(0.795, 0.152, "textos", fontsize=9, ha="left", va="center")
        area_textos = fig.add_axes([0.805, 0.115, 0.15, 0.020])
        controle_textos = Slider(area_textos, "", 0.5, 3.0,
                                  valinit=1.0, valfmt="%.2f")

        def atualizar_textos(valor):
            """Ajusta as fontes sem modificar a geometria ou os diagramas."""
            for texto, tamanho_original in textos_plotagem:
                texto.set_fontsize(valor * tamanho_original)
            ax.tick_params(axis="x", labelsize=valor * tamanho_ticks_x)
            ax.tick_params(axis="y", labelsize=valor * tamanho_ticks_y)
            fig.canvas.draw_idle()

        controle_textos.on_changed(atualizar_textos)
        sliders.append(controle_textos)

        # Exporta somente o eixo estrutural, sem os controles interativos.
        area_exportar = fig.add_axes([0.805, 0.025, 0.15, 0.04])
        botao_svg = Button(area_exportar, "Exportar SVG")

        def exportar_svg(event):
            """Salva a vista atual da análise estrutural em formato vetorial."""
            from tkinter import Tk, filedialog, messagebox
            import matplotlib as mpl

            janela = Tk()
            janela.withdraw()
            try:
                caminho = filedialog.asksaveasfilename(
                    parent=janela, title="Exportar análise estrutural",
                    defaultextension=".svg", filetypes=[("SVG", "*.svg")],
                    initialfile="analise_estrutural.svg")
                if not caminho:
                    return

                # Oculta os controles e o aviso de seleção de nós na exportação.
                artistas = [eixo for eixo in fig.axes if eixo is not ax]
                artistas += list(fig.texts)
                if not marcador_no_selecionado.get_visible():
                    artistas.append(painel_no)
                estados = [artista.get_visible() for artista in artistas]
                try:
                    for artista in artistas:
                        artista.set_visible(False)
                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    caixa = ax.get_tightbbox(renderer).transformed(
                        fig.dpi_scale_trans.inverted()).padded(0.1)
                    with mpl.rc_context({"svg.fonttype": "none"}):
                        fig.savefig(caminho, format="svg", bbox_inches=caixa)
                finally:
                    for artista, visivel in zip(artistas, estados):
                        artista.set_visible(visivel)
                    fig.canvas.draw_idle()
            except Exception as erro:
                messagebox.showerror("Exportar SVG", str(erro), parent=janela)
            finally:
                janela.destroy()

        botao_svg.on_clicked(exportar_svg)

        # Zoom ancorado no cursor (scroll) e deslocamento com botão esquerdo.
        # O pan do toolbar do Matplotlib continua disponível normalmente.
        estado_pan = {"inicio": None, "xlim": None, "ylim": None}
        def zoom_scroll(event):
            if event.inaxes is not ax or event.xdata is None or event.ydata is None:
                return
            fator = 1 / 1.2 if event.button == "up" else 1.2
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            cx, cy = event.xdata, event.ydata
            ax.set_xlim(cx + (x0-cx)*fator, cx + (x1-cx)*fator)
            ax.set_ylim(cy + (y0-cy)*fator, cy + (y1-cy)*fator)
            fig.canvas.draw_idle()

        def iniciar_pan(event):
            if event.inaxes is ax and event.button == 1 and event.xdata is not None:
                toolbar = getattr(fig.canvas, "toolbar", None)
                if toolbar is not None and getattr(toolbar, "mode", ""):
                    return
                if any(np.hypot(*(ax.transData.transform((no.x, no.y)) -
                                   np.array([event.x, event.y]))) <= 12 for no in nos):
                    return
                estado_pan.update(inicio=(event.x, event.y),
                                  xlim=ax.get_xlim(), ylim=ax.get_ylim())

        def mover_pan(event):
            if estado_pan["inicio"] is None or event.inaxes is not ax:
                return
            if event.xdata is None or event.ydata is None:
                return
            dx = (event.x - estado_pan["inicio"][0]) * (estado_pan["xlim"][1] - estado_pan["xlim"][0]) / ax.bbox.width
            dy = (event.y - estado_pan["inicio"][1]) * (estado_pan["ylim"][1] - estado_pan["ylim"][0]) / ax.bbox.height
            ax.set_xlim(*(v-dx for v in estado_pan["xlim"]))
            ax.set_ylim(*(v-dy for v in estado_pan["ylim"]))
            fig.canvas.draw_idle()

        def finalizar_pan(event):
            estado_pan["inicio"] = None

        def selecionar_no(event):
            if event.inaxes is not ax or event.button != 1:
                return
            toolbar = getattr(fig.canvas, "toolbar", None)
            if toolbar is not None and getattr(toolbar, "mode", ""):
                return
            if event.x is None or event.y is None:
                return
            candidatos = []
            for no in nos:
                px, py = ax.transData.transform((no.x, no.y))
                candidatos.append((np.hypot(px-event.x, py-event.y), no))
            distancia, no = min(candidatos, key=lambda item: item[0])
            if distancia > 12:
                return
            i = 3*(no.num-1)
            painel_no.set_text(
                f"Nó {no.num}   ux = {u[i]:.3e}   uy = {u[i+1]:.3e}   "
                f"θz = {u[i+2]:.3e} rad")
            marcador_no_selecionado.set_data([no.x], [no.y])
            marcador_no_selecionado.set_visible(True)
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_press_event", selecionar_no)
        fig.canvas.mpl_connect("scroll_event", zoom_scroll)
        fig.canvas.mpl_connect("button_press_event", iniciar_pan)
        fig.canvas.mpl_connect("motion_notify_event", mover_pan)
        fig.canvas.mpl_connect("button_release_event", finalizar_pan)

        ax.margins(0.20)
        fig._controles_analise = (checks, *sliders, botao_svg)
        if mostrar:
            plt.show()
        return fig, ax

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
                       mostrar_armadura_transversal=None,
                       decalagem_cm=None, mostrar_diagrama_decalado=None,
                       decalagem_automatica=True, ancoragem_automatica=False,
                       ancoragem_apoios=True,
                       boa_aderencia=True):
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
                mostrar_diagrama_decalado: Exibe a envoltória decaladada.
                decalagem_automatica: Estima a decalagem pelas inclinações
                    de bielas e estribos; decalagem_cm tem prioridade.
                ancoragem_apoios: Avalia ancoragem reta nas regiões dos apoios;
                    registra insuficiências geométricas sem presumir atendimento.
                ancoragem_automatica: Acrescenta comprimento necessário l_b,nec de
                    ancoragem às faixas N1/N2 após a decalagem.
                boa_aderencia: Condição de aderência adotada para as barras.
                    Deve ser confirmada para cada face e situação de concretagem.

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
        if decalagem_cm is not None:
            from math import isfinite
            if not isfinite(float(decalagem_cm)) or float(decalagem_cm) < 0:
                raise ValueError('decalagem_cm deve ser finita e não negativa.')
            self.decalagem_cm = float(decalagem_cm)
            self._decalagem_manual = True


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

        geometria_viga, m = self._desenhar_contorno(ax, nos)

        grupos = self._desenhar_modelo_estrutural(ax, barras, nos, tam_apoio, xmax, xmin, xs, ys)

        cargas = self._desenhar_cargas_distribuidas(ax, escala, grupos)

        self._desenhar_cargas_nodais(ax, escala, grupos, nos)

        solicitacoes, trechos, escala_momento, atualizar_momento = self._desenhar_momento_fletor(ax, barras, cargas, decalagem_automatica, escala, fig, grupos)

        diagramas_nv = self._desenhar_esforcos_normal_cortante(ax, escala, fig, grupos, solicitacoes, trechos)

        diagramas_armadura = self._desenhar_diagramas_armaduras(ax, barras, escala, fig, grupos, trechos)

        linhas_deformada, escala_auto, atualizar_deformada, deslocamento_no = self._desenhar_deformada_viga(ax, escala, fig, nos)

        atualizar_faixas_n1_n2, artistas_n1_n2, artistas_banzos, artistas_alertas_apoios = self._desenhar_detalhamento_longitudinal(ancoragem_apoios, ancoragem_automatica, ax, boa_aderencia, decalagem_automatica, grupos, m)

        self._desenhar_detalhamento_transversal(ax, m, grupos)

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
            'Diagrama decalado': False,
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
            'Diagrama decalado': mostrar_diagrama_decalado,
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
                'Momento fletor', 'Diagrama decalado',
                'Esforço Normal', 'Esforço cortante',
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

            # Associa cada diagrama ao eixo do seu controle de escala.
            # Os eixos são criados adiante, mas o callback já pode consultar
            # este dicionário quando o usuário alternar um checkbox.
            eixos_escalas = {}
            controles_escalonamento = []

            def atualizar_visibilidade_escalas():
                """Mostra apenas as escalas dos diagramas habilitados."""
                for nome_diagrama, eixo_escala in eixos_escalas.items():
                    eixo_escala.set_visible(bool(estados[nome_diagrama]))

            def atualizar_visibilidade_escalonamento():
                """Mostra os controles N1/N2 apenas com a armadura longitudinal."""
                visivel = bool(estados['Detalhamento longitudinal'])
                for artista in controles_escalonamento:
                    artista.set_visible(visivel)

            def alternar(rotulo):
                nome = nome_interno.get(rotulo, rotulo)
                estados[nome] = not estados[nome]
                if nome != 'Valores':
                    for artista in grupos[nome]:
                        artista.set_visible(estados[nome])
                # Ao ativar diagramas, oculta as cargas para evitar sobreposição.
                if (nome in ('Momento fletor', 'Diagrama decalado',
                             'Esforço normal',
                             'Esforço cortante', 'Armadura longitudinal',
                             'Armadura transversal') and estados[nome]
                        and estados['Cargas distribuídas']):
                    checks.set_active(nomes.index('Cargas distribuídas'))
                atualizar_valores()
                atualizar_visibilidade_escalas()
                atualizar_visibilidade_escalonamento()
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

            # ------------------------------------------------------------
            # ESCALONAMENTO N1/N2 — CONTROLES NA ÁREA DE PLOTAGEM
            # ------------------------------------------------------------
            # Os sliders ficam fixos no canto superior esquerdo dos eixos.
            # Suas coordenadas são relativas à área de plotagem, de modo que
            # o posicionamento não depende dos limites ou do zoom da viga.
            if hasattr(m, 'ap1') and hasattr(m, 'l') and hasattr(m, 'h'):
                area_n1 = ax.inset_axes(
                    [.105, .925, .225, .025],
                    transform=ax.transAxes,
                    zorder=21,
                )
                area_n2 = ax.inset_axes(
                    [.105, .865, .225, .025],
                    transform=ax.transAxes,
                    zorder=21,
                )

                for area_escalonamento in (area_n1, area_n2):
                    area_escalonamento.set_facecolor((1, 1, 1, .92))

                controle_n1 = Slider(
                    area_n1, 'N1', valmin=0, valmax=3,
                    valinit=self.niveis_n1, valstep=1, valfmt='%d',
                )
                controle_n2 = Slider(
                    area_n2, 'N2', valmin=0, valmax=3,
                    valinit=self.niveis_n2, valstep=1, valfmt='%d',
                )

                titulo_escalonamento = ax.text(
                    .015, .985, 'Escalonamento N1 / N2',
                    transform=ax.transAxes, ha='left', va='top',
                    fontsize=9, fontweight='bold', zorder=22,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.85),
                )

                def alterar_escalonamento(_):
                    """Atualiza as faixas de armadura após alterar um nível."""
                    self.niveis_n1 = int(controle_n1.val)
                    self.niveis_n2 = int(controle_n2.val)

                    atualizar_faixas_n1_n2()

                    visivel = estados['Detalhamento longitudinal']
                    for artista in (*artistas_n1_n2, *artistas_banzos,
                                    *artistas_alertas_apoios):
                        artista.set_visible(visivel)

                    fig.canvas.draw_idle()

                controle_n1.on_changed(alterar_escalonamento)
                controle_n2.on_changed(alterar_escalonamento)
                self._controles.extend((controle_n1, controle_n2))
                controles_escalonamento.extend((area_n1, area_n2, titulo_escalonamento))
                atualizar_visibilidade_escalonamento()

            if mostrar_controles_escala:
                area_escala_m = fig.add_axes([.82, .215, .14, .020])
                eixos_escalas['Momento fletor'] = area_escala_m
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
                    eixos_escalas[tipo] = area_nv
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
                    eixos_escalas[titulo] = area_arm
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

                    eixos_escalas['Deformada'] = area_slider
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

            # Respeita os estados iniciais dos checkboxes, inclusive quando
            # nenhum diagrama estiver selecionado.
            atualizar_visibilidade_escalas()

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


    def _desenhar_contorno(self, ax, nos) -> tuple:
        """
        Desenha o contorno da viga, pilares e aberturas.

        Args:

            ax: eixos Matplotlib
            nos: lista de nós

        Returns:

            geometria_viga, m.
        """

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

        return geometria_viga, m

    def _desenhar_modelo_estrutural(self, ax, barras, nos, tam_apoio, xmax, xmin, xs, ys) -> dict:
        """
        Desenha barras, nós, numeração e apoios.

        Args:

            ax: eixos Matplotlib
            barras: lista de barras
            nos: lista de nós
            tam_apoio: tam apoio
            xmax: xmax
            xmin: xmin
            xs: xs
            ys: ys

        Returns:

            grupos.
        """

        from matplotlib.patches import Circle
        from matplotlib.patches import Polygon
        from matplotlib.patches import Rectangle
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

        return grupos

    def _desenhar_cargas_distribuidas(self, ax, escala, grupos):
        """
        Desenha os carregamentos distribuídos.

        Args:

            ax: eixos Matplotlib
            escala: dimensão de referência gráfica
            grupos: dicionário de artistas gráficos

        Returns:

            cargas.
        """

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

        return cargas

    def _desenhar_cargas_nodais(self, ax, escala, grupos, nos) -> None:
        """
        Desenha forças e momentos nodais.

        Args:

            ax: eixos Matplotlib
            escala: dimensão de referência gráfica
            grupos: dicionário de artistas gráficos
            nos: lista de nós
        """

        from matplotlib.patches import Arc
        from matplotlib.patches import FancyArrowPatch
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


    def _desenhar_momento_fletor(self, ax, barras, cargas, decalagem_automatica, escala, fig, grupos) -> tuple:
        """
        Desenha o momento fletor e a envoltória decalada.

        Args:

            ax: eixos Matplotlib
            barras: lista de barras
            cargas: cargas
            decalagem_automatica: decalagem automatica
            escala: dimensão de referência gráfica
            fig: figura Matplotlib
            grupos: dicionário de artistas gráficos

        Returns:

            solicitacoes, trechos, escala_momento, atualizar_momento.
        """

        import numpy as np
        # ============================================================
        # DIAGRAMA DE MOMENTO FLETOR
        # ============================================================
        artistas_momento = []
        etiquetas_momento = []
        dados_momento = []
        extremos = []
        solicitacoes = self.modelo.portico.calcula_solicitacoes_internas(cargas)
        # Mantém os esforços das barras fictícias ocultos na apresentação.
        solicitacoes = [(b, x, np.zeros_like(N), np.zeros_like(V), np.zeros_like(M))
                        if getattr(b, "ficticia", False) else (b, x, N, V, M)
                        for b, x, N, V, M in solicitacoes]
        for barra, x_local, N, V, M in solicitacoes:
            xi, yi = barra.noi.x, barra.noi.y
            dx_b = barra.noj.x - xi
            dy_b = barra.noj.y - yi
            L = np.hypot(dx_b, dy_b)
            if L <= 0:
                continue
            c, s = dx_b/L, dy_b/L
            dados_momento.append((barra, xi, yi, c, s, x_local, M))
            if not getattr(barra, 'ficticia', False):
                extremos.extend((float(M[k]), barra.num, float(x_local[k]))
                                for k in (np.argmin(M), np.argmax(M)))

        dados_decalados = calcula_momentos_decalados(
            dados_momento, decalagem_automatica,
            self.decalagem_cm if self._decalagem_manual else None)

        Mmax = max((np.max(np.abs(d[6])) for d in dados_momento), default=0.0)
        escala_momento = 0.06 * escala / Mmax if Mmax > 1e-12 else 1.0
        objetos_momento = []
        etiquetas_extremos = []
        pontos_rotulados = set()
        extremos_ordenados = sorted(extremos, key=lambda item: item[0])
        minimo_global = extremos_ordenados[0] if extremos_ordenados else None
        maximo_global = extremos_ordenados[-1] if extremos_ordenados else None

        objetos_decalados = []
        artistas_decalados = []
        for indice, (barra, xi, yi, c, s, x_local, M) in enumerate(dados_momento):
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
            if not getattr(barra, 'ficticia', False) and abs(s) < 1e-7:
                for envoltoria in dados_decalados[indice]:
                    linha_dec, = ax.plot([], [], color='#d97706', lw=1.5,
                                         linestyle=(0, (4, 3)), zorder=9,
                                         visible=False)
                    artistas_decalados.append(linha_dec)
                    objetos_decalados.append((linha_dec, x_base, y_base,
                                              c, s, envoltoria))

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
        candidatos = [(float(M[k]), barra, float(x_local[k]))
                      for barra, _, _, _, _, x_local, M in dados_momento
                      if not getattr(barra, 'ficticia', False)
                      for k in (np.argmin(M), np.argmax(M))]
        if candidatos:
            for nome, item in (('Mmín', min(candidatos, key=lambda r: r[0])),
                               ('Mmáx', max(candidatos, key=lambda r: r[0]))):
                valor, barra, x_local = item
                c = (barra.noj.x-barra.noi.x)/barra.L
                sn = (barra.noj.y-barra.noi.y)/barra.L
                xb = barra.noi.x + c*x_local
                yb = barra.noi.y + sn*x_local
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
            for linha, xb, yb, c, s, envoltoria in objetos_decalados:
                linha.set_data(xb + s*fator*envoltoria,
                               yb - c*fator*envoltoria)
            for anot, xb, yb, c, s, valor in (
                    etiquetas_momento + etiquetas_extremos):
                anot.xy = (xb + s*fator*valor, yb - c*fator*valor)
            fig.canvas.draw_idle()

        atualizar_momento(escala_momento)
        grupos['Momento fletor'] = artistas_momento
        grupos['Diagrama decalado'] = artistas_decalados
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

        return solicitacoes, trechos, escala_momento, atualizar_momento

    def _desenhar_esforcos_normal_cortante(self, ax, escala, fig, grupos, solicitacoes, trechos) -> dict:
        """
        Desenha os diagramas de esforço normal e cortante.

        Args:

            ax: eixos Matplotlib
            escala: dimensão de referência gráfica
            fig: figura Matplotlib
            grupos: dicionário de artistas gráficos
            solicitacoes: solicitações internas das barras
            trechos: trechos geométricos das barras

        Returns:

            diagramas_nv.
        """

        import numpy as np
        # ============================================================
        # DIAGRAMAS DE ESFORCO NORMAL E CORTANTE
        # ============================================================
        diagramas_nv = {}
        for tipo, componentes, cor in (
            ('Esforço normal', (0, 3), '#7c3aed'),
            ('Esforço cortante', (1, 4), '#e07b18'),
        ):
            dados = []
            for barra, x, N, V, M in solicitacoes:
                xi, yi = barra.noi.x, barra.noi.y
                vx = barra.noj.x - xi
                vy = barra.noj.y - yi
                L = np.hypot(vx, vy)
                if L <= 0:
                    continue
                q = N if componentes == (0, 3) else V
                dados.append((barra, xi, yi, vx/L, vy/L, x, q))

            qmax = max((np.max(np.abs(d[6])) for d in dados), default=0.0)
            fator_auto = 0.06 * escala / qmax if qmax > 1e-12 else 1.0
            objetos, etiquetas, pontos_vistos = [], [], set()
            for barra, xi, yi, c, s, x, q in dados:
                xb = xi + c*x
                yb = yi + s*x
                linha, = ax.plot([], [], color=cor, lw=1.8, zorder=7,
                                 visible=False)
                preenchimento = ax.fill([], [], color=cor, alpha=0.13,
                                         zorder=3, visible=False)[0]
                objetos.append((linha, preenchimento, xb, yb, c, s, q))
            valores_nv = {id(b): (float(q[0]), float(q[-1]))
                          for b, _, _, _, _, _, q in dados}
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

        return diagramas_nv

    def _desenhar_diagramas_armaduras(self, ax, barras, escala, fig, grupos, trechos) -> dict:
        """
        Desenha os diagramas de áreas de armadura.

        Args:

            ax: eixos Matplotlib
            barras: lista de barras
            escala: dimensão de referência gráfica
            fig: figura Matplotlib
            grupos: dicionário de artistas gráficos
            trechos: trechos geométricos das barras

        Returns:

            diagramas_armadura.
        """

        import numpy as np
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

        return diagramas_armadura

    def _desenhar_deformada_viga(self, ax, escala, fig, nos) -> tuple:
        """
        Desenha a deformada e prepara sua atualização gráfica.

        Args:

            ax: eixos Matplotlib
            escala: dimensão de referência gráfica
            fig: figura Matplotlib
            nos: lista de nós

        Returns:

            linhas_deformada, escala_auto, atualizar_deformada, deslocamento_no.
        """

        import numpy as np
        deslocamento_no = None
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

            # Deslocamentos nodais para a informação de resultados.
            def deslocamento_no(no):
                i = 3 * (no.num - 1)
                return u[i:i + 3]

            max_u = np.max(np.hypot(u[0::3], u[1::3]), initial=0.0)

            if max_u > 1e-14:
                escala_auto = 0.05 * escala / max_u

            dados_deformada = portico.calcula_deformada(40)

            for barra, coordenadas, deslocamentos in dados_deformada:
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
                barra, coordenadas, deslocamentos = dados
                pontos = coordenadas + fator*deslocamentos
                linha.set_data(pontos[:, 0], pontos[:, 1])

            fig.canvas.draw_idle()

        if linhas_deformada:
            atualizar_deformada(escala_auto)



        return linhas_deformada, escala_auto, atualizar_deformada, deslocamento_no

    def _desenhar_detalhamento_longitudinal(self, ancoragem_apoios, ancoragem_automatica, ax, boa_aderencia, decalagem_automatica, grupos, m) -> tuple:
        """
        Desenha as armaduras longitudinais e prepara sua atualização.

        Args:

            ancoragem_apoios: ancoragem apoios
            ancoragem_automatica: ancoragem automatica
            ax: eixos Matplotlib
            boa_aderencia: boa aderencia
            decalagem_automatica: decalagem automatica
            grupos: dicionário de artistas gráficos
            m: modelo da viga

        Returns:

            atualizar_faixas_n1_n2, artistas_n1_n2, artistas_banzos, artistas_alertas_apoios.
        """

        atualizar_faixas_n1_n2 = None
        artistas_n1_n2, artistas_banzos, artistas_alertas_apoios = [], [], []
        # ------------------------------------------------------------------
        # DETALHAMENTO DAS ARMADURAS LONGITUDINAIS
        # ------------------------------------------------------------------
        # Detalhamento longitudinal: cores individuais e legenda fixa na área de desenho.
        grupos['Detalhamento longitudinal'] = []
        if hasattr(m, 'ap1') and hasattr(m, 'l') and hasattr(m, 'h'):
            detalhamento = Detalhamento(m)
            self.detalhamento = detalhamento
            limite_x_esquerdo = float(m.ap1)/2.0
            limite_x_direito = limite_x_esquerdo + float(m.l)
            def calcular_detalhamento():
                detalhamento.calcula_longitudinais(
                    self.niveis_n1, self.niveis_n2,
                    self.decalagem_cm if self._decalagem_manual else None,
                    decalagem_automatica, ancoragem_automatica,
                    ancoragem_apoios, boa_aderencia)
                for nome in ('armaduras_longitudinais', 'armaduras_escalonadas',
                             'armaduras_banzos', 'diagnosticos_ancoragem_apoios',
                             'diagnosticos_ancoragem_banzos'):
                    setattr(self, nome, getattr(detalhamento, nome))
            calcular_detalhamento()
            propostas = detalhamento.armaduras_longitudinais
            propostas_banzos = detalhamento.armaduras_banzos
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
            artistas_alertas_apoios = []

            def atualizar_faixas_n1_n2():
                for artista in artistas_n1_n2:
                    artista.remove()
                    if artista in grupos['Detalhamento longitudinal']:
                        grupos['Detalhamento longitudinal'].remove(artista)
                artistas_n1_n2.clear()
                for artista in artistas_alertas_apoios:
                    artista.remove()
                    if artista in grupos['Detalhamento longitudinal']:
                        grupos['Detalhamento longitudinal'].remove(artista)
                artistas_alertas_apoios.clear()
                calcular_detalhamento()
                faixas = detalhamento.armaduras_escalonadas
                novos_banzos = detalhamento.armaduras_banzos
                for artista in artistas_banzos:
                    artista.remove()
                    if artista in grupos['Detalhamento longitudinal']:
                        grupos['Detalhamento longitudinal'].remove(artista)
                artistas_banzos.clear()
                # Diagnóstico explícito: a rotina anterior apenas armazenava
                # resultados, sem apresentar avisos no terminal ou no gráfico.
                if ancoragem_automatica and ancoragem_apoios:
                    from .ancoragem_apoios import relatorio_ancoragem_apoios
                    relatorio_ancoragem_apoios(self.diagnosticos_ancoragem_apoios)
                    avisos = [d for d in self.diagnosticos_ancoragem_apoios
                              if not d.atendida_geometricamente]
                    # Marcadores e diagnósticos são recalculados a cada
                    # alteração dos controles de escalonamento.
                    import math
                    for d in avisos:
                        g = getattr(d, 'gancho_90', None)
                        if g is None:
                            continue
                        x_apoio = (limite_x_esquerdo if d.lado == 'esquerdo'
                                   else limite_x_direito)
                        y_apoio = (float(m.h) / 2.0 + 12.0
                                   if d.face == 'superior'
                                   else -float(m.h) / 2.0 - 12.0)
                        deficit_gancho = max(0.0,
                            float(d.lb_gancho_cm) - float(d.comprimento_disponivel_cm))
                        if deficit_gancho > 1e-4:
                            marcador, = ax.plot(
                                [x_apoio], [y_apoio], marker='x', markersize=10,
                                markeredgewidth=2.5, linestyle='None',
                                color='red', zorder=32, visible=False)
                            texto = ax.annotate(
                                f'! {d.posicao}: gancho — faltam {deficit_gancho:.1f} cm',
                                (x_apoio, y_apoio), xytext=(6, 7),
                                textcoords='offset points', fontsize=8,
                                color='red', weight='bold', zorder=33,
                                visible=False)
                            artistas_alertas_apoios.extend((marcador, texto))
                            grupos['Detalhamento longitudinal'].extend((marcador, texto))
                    # Alertas de insuficiência sem gancho viável.
                    for d in avisos:
                        if getattr(d, 'gancho_90', None) is not None:
                            continue
                        x_apoio = (limite_x_esquerdo if d.lado == 'esquerdo'
                                   else limite_x_direito)
                        y_apoio = (float(m.h) / 2.0 + 12.0
                                   if d.face == 'superior'
                                   else -float(m.h) / 2.0 - 12.0)
                        deficit = max(0.0, d.lb_necessario_cm - d.comprimento_disponivel_cm)
                        marcador, = ax.plot([x_apoio], [y_apoio], marker='x',
                            markersize=10, markeredgewidth=2.5,
                            linestyle='None', color='red', zorder=32, visible=False)
                        texto = ax.annotate(
                            f'! {d.posicao}: falta {deficit:.1f} cm',
                            (x_apoio, y_apoio), xytext=(6, 7),
                            textcoords='offset points', fontsize=8,
                            color='red', weight='bold', zorder=33, visible=False)
                        artistas_alertas_apoios.extend((marcador, texto))
                        grupos['Detalhamento longitudinal'].extend((marcador, texto))
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
                        f'{proposta_b.quantidade} {pos} {SIMBOLO_BITOLA}{proposta_b.bitola_mm:g} '
                        f'c = {proposta_b.comprimento_geometrico_cm:.0f} cm '
                        f'(As = {proposta_b.as_fornecida_cm2:.2f} cm²)'
                        if proposta_b is not None else '')
                # A continuidade nos furos preserva a origem N1/N2, mas
                # pode produzir barras fisicamente distintas (bitola, número,
                # comprimento ou nível). Cada tipo recebe posição própria.
                # N1/N2 ficam reservadas às barras resistentes de base.
                posicoes_ocupadas = {
                    proposta.posicao for proposta in novos_banzos
                } | {'N1', 'N2'}
                proximo_numero = 3
                identificacoes = {}
                principais = {}

                # Identifica primeiro as barras resistentes de base, evitando
                # que porta-estribos de outra bitola sejam chamados de N1/N2.
                for faixa_base in faixas:
                    if faixa_base.nivel != 0:
                        continue
                    origem = faixa_base.posicao
                    if origem not in principais:
                        principais[origem] = (
                            origem, faixa_base.face, faixa_base.quantidade,
                            round(float(faixa_base.bitola_mm), 6),
                            faixa_base.nivel,
                            round(float(faixa_base.fim_cm - faixa_base.inicio_cm), 4),
                        )

                def posicao_da_faixa(faixa):
                    """Numera tipos distintos sem modificar o aço calculado."""
                    nonlocal proximo_numero

                    chave = (
                        faixa.posicao, faixa.face, faixa.quantidade,
                        round(float(faixa.bitola_mm), 6), faixa.nivel,
                        round(float(faixa.fim_cm - faixa.inicio_cm), 4),
                    )
                    if chave in identificacoes:
                        return identificacoes[chave]

                    if chave == principais.get(faixa.posicao):
                        posicao = faixa.posicao
                    else:
                        while f'N{proximo_numero}' in posicoes_ocupadas:
                            proximo_numero += 1
                        posicao = f'N{proximo_numero}'
                        proximo_numero += 1
                        posicoes_ocupadas.add(posicao)

                    identificacoes[chave] = posicao
                    return posicao

                # Define pistas por posição efetiva, e não apenas pelo nível
                # teórico. Duas bitolas diferentes no mesmo nível precisam
                # aparecer em alturas distintas para permitir sua leitura.
                # Faixas da mesma posição compartilham a pista, inclusive
                # quando ocupam intervalos descontínuos da viga.
                pistas_por_face = {'inferior': {}, 'superior': {}}
                faixas_identificadas = []
                for faixa in faixas:
                    posicao = posicao_da_faixa(faixa)
                    faixas_identificadas.append((faixa, posicao))
                    pistas = pistas_por_face[faixa.face]
                    if posicao not in pistas:
                        pistas[posicao] = len(pistas)

                for faixa, posicao_grafica in faixas_identificadas:
                    indice_pista = pistas_por_face[faixa.face][posicao_grafica]
                    y_arm = y_faixa(faixa.face, indice_pista)
                    cor = cores_por_posicao.get(
                        faixa.posicao, '#566573'
                    )
                    # Desenho da própria barra na pista externa, incluindo
                    # a dobra e a ponta vertical, sem ganchos desconectados.
                    gancho_esq = None
                    gancho_dir = None
                    for diagnostico in self.diagnosticos_ancoragem_apoios:
                        if (diagnostico.posicao == faixa.posicao
                                and diagnostico.face == faixa.face
                                and diagnostico.gancho_90 is not None):
                            g = diagnostico.gancho_90
                            extremidade = (faixa.inicio_cm if diagnostico.lado == 'esquerdo'
                                           else faixa.fim_cm)
                            if abs(extremidade - g.x_vertical_cm) <= 6.0:
                                if diagnostico.lado == 'esquerdo':
                                    gancho_esq = g
                                else:
                                    gancho_dir = g
                    # A extremidade horizontal coincide com a tangência da
                    # curva. Nunca desenhar a linha até x_vertical, pois isso
                    # cria uma pequena ponta reta além do início do gancho.
                    x_reto_esq = (gancho_esq.x_tangencia_cm if gancho_esq is not None
                                  else faixa.inicio_cm)
                    x_reto_dir = (gancho_dir.x_tangencia_cm if gancho_dir is not None
                                  else faixa.fim_cm)

                    def pontos_gancho(g, lado):
                        """Pontos da curva do eixo da barra, da tangência à ponta."""
                        sentido_x = -1.0 if lado == 'esquerdo' else 1.0
                        angulos = [math.pi * i / 40 for i in range(21)]
                        curva_x = [g.x_tangencia_cm + sentido_x * g.raio_eixo_cm * math.sin(a)
                                   for a in angulos]
                        curva_y = [y_arm + g.sentido_vertical * g.raio_eixo_cm * (1 - math.cos(a))
                                   for a in angulos]
                        return curva_x, curva_y

                    pontos_x = [x_reto_esq, x_reto_dir]
                    pontos_y = [y_arm, y_arm]
                    if gancho_esq is not None:
                        gx, gy = pontos_gancho(gancho_esq, 'esquerdo')
                        pontos_x = [gancho_esq.x_vertical_cm, *reversed(gx), x_reto_dir]
                        pontos_y = [gy[-1] + gancho_esq.sentido_vertical * gancho_esq.ponta_vertical_cm,
                                    *reversed(gy), y_arm]
                    if gancho_dir is not None:
                        gx, gy = pontos_gancho(gancho_dir, 'direito')
                        # gx[0] é exatamente a extremidade do trecho reto.
                        pontos_x.extend(gx[1:])
                        pontos_y.extend(gy[1:])
                        pontos_x.append(gancho_dir.x_vertical_cm)
                        pontos_y.append(gy[-1] + gancho_dir.sentido_vertical * gancho_dir.ponta_vertical_cm)
                    linha, = ax.plot(pontos_x, pontos_y,
                        color=cor, lw=2.6, zorder=12, visible=False)
                    comprimento_corte = x_reto_dir - x_reto_esq
                    for g in (gancho_esq, gancho_dir):
                        if g is not None:
                            comprimento_corte += math.pi * g.raio_eixo_cm / 2 + g.ponta_vertical_cm
                    comprimento_corte = 5 * math.ceil((comprimento_corte - 1e-9) / 5)
                    rotulo = ax.annotate(
                        f'{faixa.quantidade} {posicao_grafica} {SIMBOLO_BITOLA}{faixa.bitola_mm:g} '
                        f'C={comprimento_corte:.0f} cm',
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
                texto = (f'{proposta.quantidade} {proposta.posicao} '
                         f'{SIMBOLO_BITOLA}{proposta.bitola_mm:g} '
                         f'c = {proposta.comprimento_geometrico_cm:.0f} cm '
                         f'(As = {proposta.as_fornecida_cm2:.2f} cm²)')
                anot = ax.annotate(
                    texto, xy=(0.018, 0.035), xycoords='axes fraction',
                    xytext=(0, 12 * (len(propostas_banzos) - 1 - i)),
                    textcoords='offset points', ha='left', va='bottom',
                    color=cores_por_posicao[proposta.posicao],
                    fontsize=8, family='DejaVu Sans', weight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=1.0),
                    visible=False, zorder=30, annotation_clip=False)
                grupos['Detalhamento longitudinal'].append(anot)
                anotacoes_banzos[proposta.posicao] = anot

        return atualizar_faixas_n1_n2, artistas_n1_n2, artistas_banzos, artistas_alertas_apoios

    def _desenhar_detalhamento_transversal(self, ax, m, grupos:dict) -> None:
        """
        Representa os estribos e suas legendas.

        Args:

            ax: eixos Matplotlib
            m: modelo de viga
            grupos: coleções de artistas gráficos
        """

        import numpy as np

        # ------------------------------------------------------------------
        # DETALHAMENTO DAS ARMADURAS TRANSVERSAIS
        # ------------------------------------------------------------------
        # Estribos em seção cheia e nos banzos, com relação à direita.
        grupos['Detalhamento transversal'] = []
        if hasattr(m, 'ap1') and hasattr(m, 'l') and (hasattr(m, 'furos') or hasattr(m, 'lf')):
            if not hasattr(self, 'detalhamento') or self.detalhamento.modelo is not m:
                self.detalhamento = Detalhamento(m)
            self.detalhamento.calcula_transversais()
            propostas_estribos = self.detalhamento.armaduras_transversais
            propostas_banzos_estribos = self.detalhamento.armaduras_estribos_banzos
            todos_estribos = [*propostas_estribos, *propostas_banzos_estribos]
            h = float(m.h)
            cobrimento = min(3.5, h * 0.15)
            topo, base = h / 2 - cobrimento, -h / 2 + cobrimento
            cores_estribos = ('#8050a0', '#b05b20', '#167f76', '#b03267',
                              '#4169a2', '#737c18', '#a04a36', '#6b55a8')

            # Agrupa posições com o mesmo detalhamento e a mesma altura
            # de estribo, preservando os trechos físicos de cada proposta.
            # Arredondamentos apenas estabilizam a comparação de floats.
            grupos_estribos = {}
            for proposta in todos_estribos:
                if hasattr(proposta, 'y_inferior_cm'):
                    ybaixo, yalto = (float(proposta.y_inferior_cm),
                                     float(proposta.y_superior_cm))
                else:
                    ybaixo, yalto = base, topo

                altura = yalto - ybaixo
                chave = (round(float(proposta.asw_fornecida_cm2_m), 8),
                         round(altura, 8),
                         round(float(proposta.bitola_mm), 8),
                         int(proposta.ramos),
                         round(float(proposta.espacamento_cm), 8))
                if chave not in grupos_estribos:
                    indice = len(grupos_estribos)
                    grupos_estribos[chave] = {
                        'posicao': f'E{indice + 1}',
                        'cor': cores_estribos[indice % len(cores_estribos)],
                        'quantidade': 0,
                        'proposta': proposta,
                    }
                grupo = grupos_estribos[chave]
                grupo['quantidade'] += int(proposta.quantidade)

                cor = grupo['cor']
                xs = np.linspace(proposta.inicio_cm, proposta.fim_cm,
                                 proposta.quantidade)
                for x in xs:
                    linha, = ax.plot([x, x], [ybaixo, yalto], color=cor,
                                     lw=.85, zorder=11, visible=False)
                    grupos['Detalhamento transversal'].append(linha)

                centro = (proposta.inicio_cm + proposta.fim_cm) / 2
                rotulo = ax.annotate(
                    grupo['posicao'], (centro, (ybaixo + yalto) / 2),
                    ha='center', va='center', fontsize=8, weight='bold',
                    color=cor, zorder=17, visible=False,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=1))
                grupos['Detalhamento transversal'].append(rotulo)

            # Uma linha por posição E#, somando as quantidades dos trechos.
            for i, grupo in enumerate(grupos_estribos.values()):
                proposta = grupo['proposta']
                anot = ax.annotate(
                    (f"{grupo['quantidade']} {grupo['posicao']} {SIMBOLO_BITOLA}{proposta.bitola_mm:g} "
                     f'c/{proposta.espacamento_cm:g} '
                     f'(Asws = {proposta.asw_fornecida_cm2_m:.2f} cm²/m)'),
                    xy=(.982, .035), xycoords='axes fraction',
                    xytext=(0, 12 * (len(grupos_estribos) - 1 - i)),
                    textcoords='offset points', ha='right', va='bottom',
                    color=grupo['cor'], fontsize=8, family='DejaVu Sans', weight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=1),
                    visible=False, zorder=30, annotation_clip=False)
                grupos['Detalhamento transversal'].append(anot)

            self.armaduras_transversais = propostas_estribos
            self.armaduras_estribos_banzos = propostas_banzos_estribos

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

        dados_deformada = self.modelo.portico.calcula_deformada(pontos_por_barra)

        for barra, coordenadas, deslocamentos in dados_deformada:
            pontos = coordenadas + escala_inicial*deslocamentos
            linha, = ax.plot(
                pontos[:, 0], pontos[:, 1],
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
            for linha, dados in zip(linhas, dados_deformada):
                barra, coordenadas, deslocamentos = dados
                pontos = coordenadas + fator*deslocamentos
                linha.set_data(pontos[:, 0], pontos[:, 1])

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
