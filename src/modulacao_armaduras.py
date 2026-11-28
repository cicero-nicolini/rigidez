"""Modulação dos comprimentos das armaduras."""

def comprimento_modulado(inicio_cm, fim_cm, limite_esquerdo_cm,
                          limite_direito_cm, modulo_cm=5.0):
    """Arredonda o comprimento para cima, ampliando a barra quando possível.

    Args:
        inicio_cm: Coordenada inicial da barra, em centímetros.
        fim_cm: Coordenada final da barra, em centímetros.
        limite_esquerdo_cm: Limite físico disponível à esquerda.
        limite_direito_cm: Limite físico disponível à direita.
        modulo_cm: Módulo de comprimento, em centímetros.

    Returns:
        Tupla (início, fim) com comprimento múltiplo do módulo.

    Raises:
        ValueError: Se o comprimento modular não couber nos limites.

    Notes:
        O arredondamento é sempre para cima, preservando a extensão
        necessária por decalagem e ancoragem. Não comprova a viabilidade
        física da ancoragem nas regiões de apoio ou de abertura.
    """
    from math import ceil

    inicio = float(inicio_cm)
    fim = float(fim_cm)
    minimo = float(limite_esquerdo_cm)
    maximo = float(limite_direito_cm)
    comprimento = fim - inicio
    if comprimento <= 0:
        raise ValueError('Comprimento de armadura inválido.')

    # Tolerância numérica evita arredondar 215,00000000001 para 220.
    comprimento_final = modulo_cm * ceil((comprimento - 1e-7) / modulo_cm)
    adicional = max(0.0, comprimento_final - comprimento)
    if adicional < 1e-7:
        return inicio, fim

    direita = min(adicional, max(0.0, maximo - fim))
    esquerda = adicional - direita
    if esquerda > max(0.0, inicio - minimo) + 1e-6:
        raise ValueError(
            'Não há espaço para modular o comprimento da armadura '
            f'de {comprimento:.3f} cm para {comprimento_final:.0f} cm.')

    return inicio - esquerda, fim + direita


def modular_armaduras(armaduras, limite_esquerdo_cm, limite_direito_cm):
    """Modula as coordenadas das barras e mantém seus demais atributos."""
    from dataclasses import replace

    resultado = []
    for armadura in armaduras:
        inicio, fim = comprimento_modulado(
            armadura.inicio_cm, armadura.fim_cm,
            limite_esquerdo_cm, limite_direito_cm)
        alteracoes = {'inicio_cm': inicio, 'fim_cm': fim}
        if hasattr(armadura, 'comprimento_geometrico_cm'):
            alteracoes['comprimento_geometrico_cm'] = fim - inicio
        resultado.append(replace(armadura, **alteracoes))
    return resultado


