"""Condições de contorno dos apoios de pórticos planos."""


def validar_vinculacao(classe):
    """Valida a classe de apoio (1=móvel, 2=articulado, 3=engaste)."""
    if type(classe) is not int or classe not in (1, 2, 3):
        raise ValueError('A vinculação deve ser 1 (móvel), 2 (articulado) ou 3 (engaste).')
    return classe


def aplicar_vinculacao(no, classe):
    """Configura as restrições Tx, Ty e Rz de um nó de pórtico plano."""
    validar_vinculacao(classe)
    no.Tx = classe >= 2
    no.Ty = True
    no.Rz = classe == 3
