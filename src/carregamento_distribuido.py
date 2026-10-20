from barra import Barra
import numpy as np

class Carregamento_distribuido:
    """Representa as propriedades e parâmetros referente ao objeto carregamento distribuído."""
    
    def __init__(self, a:float, lw:float, w1:float, w2:float, barra:Barra) -> None:
        """
        Instancia um objeto da classe Carregamento_distribuido


        Args:

            a: distância do nó inicial da barra até o inicio do carregamento
            lw: comprimento do carregamento
            w1: modulo inicial do carregamento
            w2: modulo final do carregamento
            barra: objeto barra
        """

        self.a = a
        self.lw = lw
        self.w1 = w1
        self.w2 = w2
        self.barra = barra
        
    def calcula_fepl(self) -> None:
        """Calcula o vetor de forças de engastamento perfeito no sistema local da barra."""

        L = self.barra.L
        b = L - self.lw - self.a
        wm = (self.w1 + self.w2)/2.0
        wd = self.w2 - self.w1
        s1 = 10.0*((L*L + self.a*self.a)*(L+self.a)-(self.a*self.a+b*b)*(self.a-b)-(L*b)*(L+b)-self.a**3.0)
        s2 = self.lw*(L*(2.0*L+self.a+b)-3.0*(self.a-b)*(self.a-b)-2.0*self.a*b)
        s3 = 120.0*self.a*b*(self.a + self.lw)+10.0*self.lw*(6.0*self.a*self.a+4.0*L*self.lw-3.0*self.lw*self.lw)
        s4 = 10.0*L*self.lw*self.lw-10.0*self.lw*self.a*(L-3.0*b)-9.0*self.lw**3

        rb = (self.lw*(s1 * wm + s2 * wd ))/(20.0*L**3)
        ra = self.lw*wm-rb
        mb = -(self.lw*(s3*wm+s4*wd))/(120.0*L*L)
        ma = -mb+rb*L-self.a*self.lw*wm-(self.lw*self.lw*(2.0*self.w2+self.w1))/6.0
        
        self.barra.fepl[1] += ra
        self.barra.fepl[2] += ma
        self.barra.fepl[4] += rb
        self.barra.fepl[5] += mb

    def calcula_fep(self) -> None:
        """Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra."""        

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl
