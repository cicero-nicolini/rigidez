from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Exemplo 1")
print("#########################################################################################################")

# Definição dos nós
no1 = No(1,0.0,75.0)
no2 = No(2,100.0,75.0)
no3 = No(3,200.0,0.0)
nos = [no1, no2, no3]

# Aplicação das cargas nodais
no2.Fy = -10.0
no2.Mz = -1000.0

# Aplicação das restrições nodais
no1.Tx = True
no1.Ty = True
no1.Rz = True
no3.Tx = True
no3.Ty = True
no3.Rz = True

# Definição das barras e propriedades
barra1 = Barra(1.0,no1,no2,10000.0,10.0,1000.0)
barra2 = Barra(2.0,no2,no3,10000.0,10.0,1000.0)
barras = [barra1, barra2]

# Definição dos carregamentos nas barras
carregamento1 = Carregamento_distribuido(0,100,0.24,0.24,barra1)
carregamento2 = Carregamento_pontual(62.5,12.0,16.0,barra2)

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

# Verificação Fep
print("fep barra1")
print_vetor(barra1.fep)
print("fep barra2")
print_vetor(barra2.fep)

# Monta vetor de cargas nodais global
portico.monta_fnos()
print("fnos")
print_vetor(portico.fnos)

# Aplica condições de contorno no vetor de cargas nodais global
portico.aplica_cc_fnos()
print("f01")
print_vetor(portico.fnos)

# Calcula vetor de deslocamento da estrutura
portico.calcula_deslocamentos()
print("deslocamentos")
print_vetor(portico.u)
portico.calcula_solicitacoes_internas_nodais()

# Verificação vetor f
print("f barra1")
print_vetor(barra1.f)
print("f barra2")
print_vetor(barra2.f)

# Calcula solicitacoes internas nas barras
portico.calcula_solicitacoes_internas_nodais()
print("solicitacoes barra1")
print_vetor(barra1.fl)
print("solicitacoes barra2")
print_vetor(barra2.fl)

# Calcula vetor de reações da estrutura
portico.calcula_reacoes()
print("reacoes")
print_vetor(portico.R)
