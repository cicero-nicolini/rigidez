
class No:
    """Representa as propriedades e parâmetros referente ao objeto nó."""

    def __init__(self, num:int, x:float, y:float) -> None:
        """
        Instancia um objeto da classe Nó


        Args:

            num: número do nó
            x: coordenada x do nó
            y: coordenada y do nó
        """

        self.num = int(num) 
        self.x = x
        self.y = y
        self.Fx = 0.0
        self.Fy = 0.0
        self.Mz = 0.0
        self.Tx = False
        self.Ty = False
        self.Rz = False