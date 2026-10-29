"""Definição geométrica de uma abertura retangular em uma viga."""
from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class Furo:
    """Abertura retangular; dimensões na unidade adotada pelo modelo (cm).

    lf: distância do início do vão (face interna do pilar esquerdo)
        até o centro longitudinal da abertura.
    b: comprimento horizontal da abertura.
    a: altura livre da abertura.
    c: distância entre a face inferior da viga e a abertura.
    """
    lf: float
    b: float
    a: float
    c: float

    def __post_init__(self):
        for nome in ('lf', 'b', 'a', 'c'):
            valor = getattr(self, nome)
            if not isfinite(valor):
                raise ValueError(f'{nome} deve ser finito.')
        if self.b <= 0 or self.a <= 0 or self.c <= 0:
            raise ValueError('b, a e c devem ser positivos.')

    @property
    def inicio(self) -> float:
        return self.lf - self.b / 2.0

    @property
    def fim(self) -> float:
        return self.lf + self.b / 2.0

    def validar_na_viga(self, l: float, h: float, bw: float = 0.0) -> None:
        """Valida abertura e zona de transição com largura bw/2 em cada lado."""
        if self.inicio - bw / 2 <= 0 or self.fim + bw / 2 >= l:
            raise ValueError('Abertura ou sua zona de transição excede o vão da viga.')
        if self.a + self.c >= h:
            raise ValueError('A abertura deve deixar um banzo superior positivo.')

    def sobrepoe(self, outro: 'Furo', folga: float = 0.0) -> bool:
        """Verifica sobreposição horizontal considerando folga mínima."""
        return self.inicio < outro.fim + folga and outro.inicio < self.fim + folga
