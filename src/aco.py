class Aco:
    """Representa as propriedades e parâmetros referente ao aço."""

    def __init__(self, fyk:float, fywk:float, γs:float = 1.15, Es:float = 21000.0) -> None:
        """
        Instancia um objeto da classe Aco


        Args:

            fyk: resistência característica do aço [MPa]
            γs: coeficiente de segurança do aço
            Es: módulo de elasticidade do aço
        """

        self.fyk = fyk
        self.fywk = fywk 
        self.γs = γs
        self.Es = Es #Unidade compatibilizada para testes de MPa para kN/cm²
    
    @property
    def fyd(self) -> float:
        fyd = (self.fyk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fyd

    @property
    def fywd(self) -> float:
        fywd = (self.fywk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fywd
    
    @property
    def εyd(self) -> float:
        fyd = (self.fyk/self.γs)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        εyd = fyd/self.Es
        return εyd