
class Secao:
    """Representa as propriedades e parâmetros referente a seção transversal."""

    def __init__(self,nome:str, b:float, h:float, E:float = None, d:float = None, d_linha:float = None) -> None:
        """
        Instancia um objeto da classe Secao


        Args:

            nome: nome da seção
            b: largura da seção
            h: altura da seção
            E: módulo de elasticidade da seção
            d: altura útil da armadura tracionada
            d_linha: altura útil da armadura comprimida 
        """

        self.nome = nome
        self.b = b
        self.h = h
        self.d = d
        self.d_linha = d_linha
        self.E = E
        self.A = b*h
        self.I = b*h**3/12
        self.yt = h/2 #distância do centro de gravidade até a fibra mais tracionada
        self.W0 = (b*h**3/12)/(h/2) #módulo de resistência à flexão da seção bruta de concreto 
                                    #em relação à fibra mais tracionada
