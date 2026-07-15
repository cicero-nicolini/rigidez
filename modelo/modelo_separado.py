from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Exemplo Modelo")
print("#########################################################################################################")


# Definição dos nós
no1 = No(1,0.0,0.0)
no2 = No(2,10.0,0.0)
no3 = No(3,90.0,0.0)
no4 = No(4,90.0,-20.0)
no5 = No(5,90.0,20.0)
no6 = No(6,130.0,0.0)
no7 = No(7,130.0,-20.0)
no8 = No(8,130.0,20.0)
no9 = No(9,230.0,0.0)
no10 = No(10,240.0,0.0)
nos = [no1, no2, no3, no4, no5, no6, no7, no8, no9, no10]

# Aplicação das restrições nodais
no1.Tx = True
no1.Ty = True
no1.Rz = True
no10.Tx = True
no10.Ty = True
no10.Rz = True

# Definição das barras e propriedades
barra1 = Barra(1.0,no1,no2,280000,1000.0,2.0833e5)
barra2 = Barra(2.0,no2,no3,280000,1000.0,2.0833e5)
barra3 = Barra(3.0,no4,no7,280000,200.0,1666.6667)
barra4 = Barra(4.0,no5,no8,280000,200.0,1666.6667)
barra5 = Barra(5.0,no6,no9,280000,1000.0,2.0833e5)
barra6 = Barra(6.0,no9,no10,280000,1000.0,2.0833e5)
barra7 = Barra(7.0,no3,no4,280000,1000000.0,2.0833e8)
barra8 = Barra(8.0,no3,no5,280000,1000000.0,2.0833e8)
barra9 = Barra(9.0,no6,no7,280000,1000000.0,2.0833e8)
barra10 = Barra(10.0,no6,no8,280000,1000000.0,2.0833e8)
barras = [barra1, barra2, barra3, barra4, barra5, barra6, barra7, barra8, barra9, barra10]

# Definição dos carregamentos nas barras
carregamento1 = Carregamento_distribuido(0,80.0,10.0,10.0,barra2)
carregamento2 = Carregamento_distribuido(0,40.0,10.0,10.0,barra4)
carregamento3 = Carregamento_distribuido(0,100.0,10.0,10.0,barra5)


# Definição da estrutura
portico = Estrutura(nos,barras)

# Monta matriz de rigidez global
portico.monta_k()
print("k")
print_matriz(portico.k)

# Monta matriz de rigidez global com as condições de contorno
portico.monta_k01()
print("k01")
print_matriz(portico.k01)

# Calcula as forças de engastamento perfeito
carregamento1.calcula_fepl()
carregamento1.calcula_fep()
carregamento2.calcula_fepl()
carregamento2.calcula_fep()
carregamento3.calcula_fepl()
carregamento3.calcula_fep()

# Monta vetor de cargas nodais global
portico.monta_fnos()
print("fnos")
print_vetor(portico.fnos)

# Aplica condições de contorno no vetor de cargas nodais global
print("f01")
portico.aplica_cc_fnos()
print_vetor(portico.fnos)

# Calcula vetor de deslocamento da estrutura
print("deslocamentos")
portico.calcula_deslocamentos()
print_vetor(portico.u)

# Calcula solicitacoes internas nas barras
portico.calcula_solicitacoes_internas_nodais()
print("solicitacoes barra1")
print_vetor(barra1.fl)
print("solicitacoes barra2")
print_vetor(barra2.fl)
print("solicitacoes barra3")
print_vetor(barra3.fl)
print("solicitacoes barra4")
print_vetor(barra4.fl)
print("solicitacoes barra5")
print_vetor(barra5.fl)
print("solicitacoes barra6")
print_vetor(barra6.fl)
print("solicitacoes barra7")
print_vetor(barra7.fl)
print("solicitacoes barra8")
print_vetor(barra8.fl)
print("solicitacoes barra9")
print_vetor(barra9.fl)
print("solicitacoes barra10")
print_vetor(barra10.fl)

# Calcula vetor de reações da estrutura
print("reacoes")
portico.calcula_reacoes()
print_vetor(portico.R)