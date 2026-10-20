from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Modelo")
print("#########################################################################################################")

#modelo_teste1 = Modelo(20.0, 50.0, 220.0, 100.0, 30.0, 40.0, 10.0, 20.0, 20.0, 10.0, 25.0, 500.0, 500.0, 1)
modelo = Modelo(bw = 30.0,
                h = 80.0,
                l = 940.0,
                lf = 150.0+60.0,
                a = 25.0,
                b = 120.0,
                c = 30.0,
                ap1 = 40.0,
                ap2 = 40.0,
                w = 0.45+0.06+0.02-0.0023,
                fck = 35.0,
                fyk = 500.0,
                fywk = 500.0,
                modelo = 1,
                α = 90,
                θ = 30,
                αE = 1)

modelo.desenha_area_armaduras()
