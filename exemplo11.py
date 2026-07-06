from rigidez import *

print("#########################################################################################################")
print("Teste Exemplo 11")
print("#########################################################################################################")


# Definição dos nós
no1 = No(1,0.0,0.0)
no2 = No(2,100.0,0.0)
no3 = No(3,200.0,50.0)
no4 = No(4,0.0,100.0)
no5 = No(5,100.0,100.0)
no6 = No(6,200.0,100.0)
nos = [no1, no2, no3, no4, no5, no6]

# Aplicação das cargas nodais
no4.Fy = -2
no4.Mz = -5
no6.Fy = -2
no6.Mz = -5

# Aplicação das restrições nodais
no1.Tx = True
no1.Ty = True
no1.Rz = True
no2.Tx = True
no2.Ty = True
no2.Rz = True
no3.Tx = True
no3.Ty = True
no3.Rz = True

# Definição das barras e propriedades
barra1 = Barra(1.0,no1,no4,2.5e2,200.0,6670.0)
barra2 = Barra(2.0,no2,no5,2.5e2,200.0,6670.0)
barra3 = Barra(3.0,no3,no6,2.5e2,200.0,6670.0)
barra4 = Barra(4.0,no4,no5,2.5e2,200.0,6670.0)
barra5 = Barra(5.0,no5,no6,2.5e2,200.0,6670.0)
barras = [barra1, barra2, barra3, barra4, barra5]

# Definição dos carregamentos nas barras
carregamento1 = Carregamento_distribuido(0,100,0.1,0.1,barra4)

# Definição da estrutura
portico = Estrutura(nos,barras)

# Monta matriz de rigidez global
print("k")
portico.monta_k()
print(portico.k)

# Monta matriz de rigidez global com as condições de contorno
portico.monta_k01()

# Calcula as forças de engastamento perfeito
carregamento1.calcula_fepl()
carregamento1.calcula_fep()

# Monta vetor de cargas nodais global
portico.monta_fnos()

# Aplica condições de contorno no vetor de cargas nodais global
print("f01")
portico.aplica_cc_fnos()
print(portico.fnos)

# Calcula vetor de deslocamento da estrutura
print("deslocamentos")
portico.calcula_deslocamentos()
print(portico.u)

portico.calcula_solicitacoes_internas_nodais()
print("solicitacoes")
print(barra1.fl)
print(barra2.fl)
print(barra3.fl)
print(barra4.fl)
print(barra5.fl)

# Calcula vetor de reações da estrutura
print("reacoes")
portico.calcula_reacoes()
print(portico.R)