"""Funções auxiliares de apresentação numérica."""

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


