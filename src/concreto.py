import math

class Concreto:
    """Representa as propriedades e parâmetros referente ao concreto."""

    def __init__(self, fck: float, γc: float = 1.4, αE: float = 1.0) -> None:
        """
        Instancia um objeto da classe Concreto


        Args:

            fck: resistência característica do concreto [MPa]
            γc: coeficiente de segurança do concreto
            αE: coeficiente de ajuste do módulo de elasticidade

            Valores de αE:

            αE = 1.2 para basalto e diabásio
            αE = 1.0 para granito e gnaisse
            αE = 0.9 para calcário
            αE = 0.7 para arenito

        Raises:
            ValueError: Se o fck estiver fora do intervalo de 20 a 90 MPa.
            ValueError: Caso αE não corresponda aos valores indicados.
        """

        if fck > 90 or fck < 20:
            raise ValueError("fck deve ser entre 20 e 90")
        
        if αE not in [0.7, 0.9, 1.0, 1.2]:
            raise ValueError("αE deve ser 0.7, 0.9, 1.0 ou 1.2")
        
        self.fck = fck
        self.γc = γc
        self.αE = αE
    
    @property
    def fct_m(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fct_m

    @property
    def fctk_sup(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_sup = 1.3*fct_m
        return fctk_sup
        
    @property
    def fctk_inf(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_inf = 0.7*fct_m
        return fctk_inf

    @property
    def fctd(self) -> float:
        fct_m = (0.3*self.fck**(2/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            fct_m = (2.12*math.log(1 + 0.1*(self.fck + 8)))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        fctk_inf = 0.7*fct_m
        fctd = fctk_inf/self.γc
        return fctd

    @property
    def αv2(self) -> float:
        αv2 = (1 - (self.fck/250))
        return αv2 

    @property
    def fcd(self) -> float:
        fcd = (self.fck/self.γc)/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return fcd
    
    @property
    def ηc(self) -> float:
        ηc = 1.0
        if 40 < self.fck <= 90:
            ηc = (40/self.fck)**(1/3)
        return ηc
    
    @property
    def Eci(self) -> float:
        Eci = (5600*self.αE*math.sqrt(self.fck))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        if 50 < self.fck <= 90:
            Eci = (21.5*10**3*self.αE*((self.fck/10)+1.25)**(1/3))/10 #Compatibilização de unidades para testes de MPa para kN/cm²
        return Eci

    @property
    def Ecs(self) -> float:
        alphai = 0.8+0.2*self.fck/80.0
        Ecs = alphai*self.Eci
        return Ecs

    @property
    def αc(self) -> float:
        αc = 0.85
        if 50 < self.fck <= 90:
            αc = 0.85*(1-(self.fck-50)/200)
        return αc

    @property
    def λ(self) -> float:
        λ = 0.8
        if 50 < self.fck <= 90:
            λ = 0.8 - ((self.fck-50)/400)
        return λ

    @property
    def εc2(self) -> float:
        εc2 = 0.002
        if 50 < self.fck <= 90:
            εc2 = 0.002 + 0.000085*((self.fck-50)**0.53)
        return εc2

    @property
    def εcu(self) -> float:
        εcu = 0.0035
        if 50 < self.fck <= 90:
            εcu = 0.0026 + 0.035*((90-self.fck)/100)**4
        return εcu