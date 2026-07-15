from pathlib import Path
import sys

RAIZ_PROJETO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ_PROJETO))

from rigidez import *

print("#########################################################################################################")
print("Teste Exemplo 1 Modelo")
print("#########################################################################################################")


# Definição dos nós
no1 = No(1,0.0,0.0)
no2 = No(2,80.0,0.0)
nos = [no1, no2]

# Aplicação das restrições nodais
no1.Tx = True
no1.Ty = True
no1.Rz = True
no2.Tx = True
no2.Ty = True
no2.Rz = True

# Definição das barras e propriedades
barra1 = Barra(1.0,no1,no2,280000,1000.0,2.0833e5)
barras = [barra1]

# Definição dos carregamentos nas barras
carregamento1 = Carregamento_distribuido(0,80.0,10.0,10.0,barra1)

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

# Verificação Fep
print("fep barra1")
print_vetor(barra1.fep)

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

# Calcula vetor de reações da estrutura
print("reacoes")
portico.calcula_reacoes()
print_vetor(portico.R)