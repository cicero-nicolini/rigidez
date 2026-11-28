                                                         
import numpy as np

from src import barra

class Estrutura:
    """Representa as propriedades e parâmetros referente ao objeto estrutura."""

    def __init__(self, nos:list, barras:list) -> None:
        """
        Instancia um objeto da classe Estrutura


        Args:

            nos: lista de objetos nós
            barras: lista de objetos barras
        """

        self.nos = nos
        self.barras = barras
        self.nnos = len(nos)
        self.fnos = np.zeros((self.nnos*3))
        self.f = np.zeros((self.nnos*3))
        self.u = np.zeros((self.nnos*3))
        self.k = np.zeros((self.nnos*3,self.nnos*3))
        self.k01 = np.zeros((self.nnos*3,self.nnos*3))
        self.R = np.zeros((self.nnos*3))
              
    def monta_k(self) -> None:
        """Monta a matriz de rigidez da estrutura."""        
        
        for barra in self.barras:
            barra.comprimento_barra()
            barra.calcula_klocal()
            barra.calcula_r()
            barra.calcula_k()
            barra.monta_q()
            for j in range(0,6):
                for jk in range(0,6):
                    self.k[barra.q[j], barra.q[jk]] += barra.k[j,jk]

    def monta_k01(self) -> None:
        """Aplica as condições de contorno na matriz da estrutura para calcular os deslocamentos."""     
    
        self.k01 = self.k
        for no in self.nos:
            if no.Tx:
                gdl = 3*no.num-3
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0
            if no.Ty:
                gdl = 3*no.num-2
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0               
            if no.Rz:
                gdl = 3*no.num-1
                self.k01[gdl, :] = 0.0
                self.k01[:,gdl] = 0.0
                self.k01[gdl,gdl] = 1.0           
        
    def monta_fnos(self) -> None:
        """Monta vetor de forças necessário para calcular deslocamentos na estrutura."""    

        for no in self.nos:
            self.fnos[3*no.num-3] = no.Fx
            self.fnos[3*no.num-2] = no.Fy
            self.fnos[3*no.num-1] = no.Mz 
            
        for barra in self.barras:
            for j in range(0,6):
                self.fnos[barra.q[j]] += -barra.fep[j]
                
    def aplica_cc_fnos(self) -> None:
        """Aplica as condições de contorno no vetor fnos para calcular deslocamentos na estrutura."""      
    
        for no in self.nos:
            if no.Tx:
                self.fnos[3*no.num-3] = 0
            if no.Ty:
                self.fnos[3*no.num-2] = 0
            if no.Rz:
                self.fnos[3*no.num-1] = 0

    def calcula_deslocamentos(self) -> None:
        """Calcula o vetor de deslocamentos da estrutura."""     
        
        self.u = np.linalg.inv(self.k01)@ self.fnos

    def calcula_deformada(self, pontos_por_barra:int = 40) -> list:
        """
        Calcula os deslocamentos ao longo das barras da estrutura.

        Args:

            pontos_por_barra: número de pontos de interpolação por barra

        Returns:

            Lista de tuplas com a barra, as coordenadas originais e os
            deslocamentos interpolados no sistema global (matrizes n x 2).

        Raises:

            ValueError: número de pontos inferior a dois, vetor de
                deslocamentos inválido ou barra com comprimento nulo.
        """

        if pontos_por_barra < 2:
            raise ValueError("Número de pontos inferior a dois.")

        u = np.asarray(self.u, dtype=float).ravel()
        if u.size != 3*self.nnos or not np.all(np.isfinite(u)):
            raise ValueError("Vetor de deslocamentos inválido ou incompatível.")

        deformada = []

        # Coordenada adimensional ao longo da barra
        r = np.linspace(0.0, 1.0, pontos_por_barra)

        for barra in self.barras:
            barra.comprimento_barra()
            L = barra.L
            if L <= 0.0:
                raise ValueError("Barra com comprimento nulo.")

            barra.monta_q()
            barra.calcula_r()

            # Transforma os deslocamentos nodais para o sistema local
            ul = barra.r @ u[barra.q]
            u1, v1, t1, u2, v2, t2 = ul

            # Interpolação linear dos deslocamentos axiais
            ua = (1-r)*u1 + r*u2

            # Funções de interpolação cúbica de Hermite
            H1 = 1 - 3*r**2 + 2*r**3
            H2 = L*(r - 2*r**2 + r**3)
            H3 = 3*r**2 - 2*r**3
            H4 = L*(-r**2 + r**3)

            # Interpolação dos deslocamentos transversais
            va = H1*v1 + H2*t1 + H3*v2 + H4*t2

            # Coordenadas originais no sistema global
            R = barra.r[:2,:2]
            pi = np.array([barra.noi.x, barra.noi.y])
            coordenadas = pi + np.column_stack((L*r, np.zeros_like(r))) @ R

            # Transforma os deslocamentos interpolados para o sistema global
            deslocamentos = np.column_stack((ua, va)) @ R

            deformada.append((barra, coordenadas, deslocamentos))

        return deformada

    def calcula_solicitacoes_internas_nodais(self) -> None:
        """Calcula as forças no sistema local para as barras = solicitações."""
        
        for barra in self.barras:
            barra.monta_q()
            barra.calcula_r()

            # Monta o vetor u para as barras
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]

            #Calcula o vetor f para as barras

            barra.f = barra.k @ barra.u + barra.fep

            #Transforma o vetor f para o sistema local para as barras = solicitacoes
                 
            barra.fl = barra.r @ barra.f 

    def calcula_ASL_barra(self) -> None:
        """Calcula as armaduras considerando as solicitações nos extremos da barra."""
        for barra in self.barras:

            # Barras fictícias não são dimensionadas.
            if barra.secao.d <= 0.0:
                barra.ASL = np.zeros(2)
                barra.ASL2 = np.zeros(2)
                barra.ASW_S = np.zeros(2)
                continue

            barra.calcula_vetor_ASL()
            barra.calcula_vetor_ASW_S()

    def calcula_reacoes(self) -> None:
        """Calcula um vetor com as reações da estrutura."""      
        
        for barra in self.barras:
            barra.monta_q()

            #Monta o vetor u para as barras
            
            for i in range(0,6):
                barra.u[i] = self.u[barra.q[i]]
            
            #Calcula o vetor f para as barras       
            
            barra.f = barra.k @ barra.u + barra.fep

            if barra.noi.Tx:
                self.R[3*barra.noi.num-3] += barra.f[0]
            if barra.noi.Ty:
                self.R[3*barra.noi.num-2] += barra.f[1]
            if barra.noi.Rz:
                self.R[3*barra.noi.num-1] += barra.f[2]

            if barra.noj.Tx:
                self.R[3*barra.noj.num-3] += barra.f[3]
            if barra.noj.Ty:
                self.R[3*barra.noj.num-2] += barra.f[4]
            if barra.noj.Rz:
                self.R[3*barra.noj.num-1] += barra.f[5]

        for no in self.nos:

            if no.Tx:
                self.R[3*no.num-3] += -no.Fx
            if no.Ty:
                self.R[3*no.num-2] += -no.Fy
            if no.Rz:
                self.R[3*no.num-1] += -no.Mz
