                                                         
import numpy as np

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

            if barra.num == 7 or barra.num == 8 or barra.num == 9 or barra.num == 10:
                pass
            else:
                barra.calcula_vetor_ASL()

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
