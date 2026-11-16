"""Pós-processamento de esforços internos de barras de pórtico plano."""

import numpy as np


def calcular_diagramas_barra(barra, carregamentos=(), n_pontos=201):
    """Reconstrói N, V e M por equilíbrio no sistema local da barra.

    Args:
        barra: Barra com esforços nodais locais ``fl`` já calculados.
        carregamentos: Cargas pontuais ou distribuídas associadas à barra.
        n_pontos: Número mínimo de pontos da discretização.

    Returns:
        Dicionário com vetores ``x``, ``N``, ``V`` e ``M``.
        Em forças pontuais, há coordenadas repetidas para mostrar os saltos.

    Raises:
        ValueError: Para dados geométricos ou cargas incompatíveis.
    """
    L = float(barra.L)
    if L <= 0 or n_pontos < 2:
        raise ValueError('Comprimento e número de pontos devem ser positivos.')
    fl = np.asarray(barra.fl, dtype=float)
    if fl.shape != (6,) or not np.all(np.isfinite(fl)):
        raise ValueError('O vetor barra.fl deve conter seis valores finitos.')

    cargas = [c for c in carregamentos if c.barra is barra]
    pontos = [0.0, L]
    for c in cargas:
        if hasattr(c, 'lw'):
            a, b = float(c.a), float(c.a + c.lw)
            if not (0 <= a < b <= L):
                raise ValueError('Intervalo de carga distribuída fora da barra.')
            pontos.extend([a, b])
        else:
            a = float(c.a)
            if not 0 <= a <= L:
                raise ValueError('Força concentrada fora da barra.')
            pontos.append(a)
    xs = sorted(set(np.linspace(0, L, int(n_pontos)).tolist() + pontos))

    def esforcos(x, incluir_pontuais):
        # Esforços internos no lado imediatamente à direita do nó inicial.
        N, V, M = -fl[0], -fl[1], -fl[2]
        for c in cargas:
            a = float(c.a)
            if hasattr(c, 'lw'):
                d = min(max(x-a, 0.0), float(c.lw))
                if d <= 0:
                    continue
                k = (float(c.w2)-float(c.w1))/float(c.lw)
                # Carga distribuída positiva no cadastro atua em -y local.
                # Integrais de q e q*(x-s), no trecho efetivamente carregado.
                resultante = float(c.w1)*d + k*d*d/2
                primeiro = (x-a)*resultante - (float(c.w1)*d*d/2 + k*d**3/3)
                V += resultante
                M -= primeiro
            elif x > a or (incluir_pontuais and x >= a):
                # Px positivo no eixo x local; Py positivo para baixo.
                N -= float(c.Px)
                V += float(c.Py)
                M -= float(c.Py)*(x-a)
        # Integração da parcela sem carregamentos.
        M -= (-fl[1])*x
        return N, V, M

    # Acrescenta as posições de V(x)=0, que são candidatas a extremos
    # interiores do momento. Não atravessa descontinuidades de força pontual.
    pontos_criticos = []
    for esquerda, direita in zip(xs[:-1], xs[1:]):
        ve = esforcos(esquerda, True)[1]
        vd = esforcos(direita, False)[1]
        if not np.isfinite(ve) or not np.isfinite(vd):
            continue
        if ve * vd >= 0 or np.isclose(ve, 0) or np.isclose(vd, 0):
            continue
        a_raiz, b_raiz = esquerda, direita
        for _ in range(48):
            meio = (a_raiz + b_raiz) / 2
            vm = esforcos(meio, True)[1]
            if ve * vm > 0:
                a_raiz, ve = meio, vm
            else:
                b_raiz = meio
        pontos_criticos.append((a_raiz + b_raiz) / 2)
    xs = sorted(set(xs + pontos_criticos))

    x_out, valores = [], []
    for x in xs:
        concentrada = any(not hasattr(c, 'lw') and np.isclose(x, c.a, atol=1e-10) for c in cargas)
        if concentrada:
            x_out.append(x)
            valores.append(esforcos(x, False))
        x_out.append(x)
        valores.append(esforcos(x, True))
    matriz = np.asarray(valores)
    return {'x': np.asarray(x_out), 'N': matriz[:, 0], 'V': matriz[:, 1], 'M': matriz[:, 2]}
