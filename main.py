from rigidez import Estrutura, No, Barra, CargaDistribuida, CargaPontual

meu_nos = [
    No(60,135,True),
    No(160,135),
    No(260,60,True),
    #No(20,20),
]

c1 = CargaDistribuida(0, 100, 0.24, 0.24)
c2 = CargaPontual(62.5, 12, 16)

meu_barras = [
    Barra(0,1,[c1]),
    Barra(1,2,[c2]),
    #Barra(2,3),
    #Barra(0,2),
]

e = Estrutura(meu_nos, meu_barras)
matriz_estrutura, forcas = e.matriz_estrutura()
