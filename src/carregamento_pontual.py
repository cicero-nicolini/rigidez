from .barra import Barra
import numpy as np

class Carregamento_pontual:
    """Representa as propriedades e parâmetros referente ao objeto carregamento pontual."""
   
    def __init__(self, a:float, Px:float, Py:float, barra:Barra) -> None:
        """
        Instancia um objeto da classe Carregamento_pontual


        Args:

            a: distância do nó inicial da barra até o ponto de aplicação da força
            Px: componente x da força
            Py: componente y da força
            barra: objeto barra
        """
        
        self.a = a
        self.Px = Px
        self.Py = Py
        self.barra = barra

    def calcula_fepl(self) -> None:
        """Calcula o vetor de forças de engastamento perfeito no sistema local da barra."""
        
        L = self.barra.L
        b = L - self.a
        sa = L + 2*self.a
        sb = L + 2*b

        ha = -(self.Px*b*b*sa)/(L**3)
        hb = -(self.Px*self.a*self.a*sb)/(L**3)
        ra = (self.Py*b*b*sa)/(L**3)
        rb = (self.Py*self.a*self.a*sb)/(L**3)
        ma = (self.Py*self.a*b*b)/(L*L)
        mb = -(self.Py*self.a*self.a*b)/(L*L)

        self.barra.fepl[0] += ha
        self.barra.fepl[1] += ra
        self.barra.fepl[2] += ma
        self.barra.fepl[3] += hb
        self.barra.fepl[4] += rb
        self.barra.fepl[5] += mb

    def calcula_fep(self) -> None:
        """Transforma o vetor de forças de engastamento perfeito no sistema local para o global na barra."""

        self.barra.fep = np.linalg.inv(self.barra.r)@ self.barra.fepl
