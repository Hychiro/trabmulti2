import numpy as np
import math
import matplotlib.pyplot as plt

class PorositiyProblem():
    def __init__(self):
        self.krw0__ = 0.7737 # ajustavel
        self.krg0__ = 0.20 # ajustavel
        self.lamb__ = 0.52 # ajustavel
        self.swc__ = 0.2 # ajustavel
        self.sw__ = 1.0 # vai variar na conta, mas é o valor inicial da saturação
        self.sg__ = 0.0 # constante, sg é o resultado do sw-swc se me lembro bem (saturação de liquido - saturação conata = saturação do gas) já q sgr pode ser considerado 0
        self.mrf__ = 1.0 # valor referente ao sufactante 
        self.phi__ = 0.3 # valor da porosidade
        self.volW__ = 159.893 # volume de agua dentro do experimento (resultado da massaFinal - massaInicial do experimento)
        self.volt__ = 319.785 # volume total do ambiente
        
        #Condicao inicial
        self.Sw_0 = 1.0 #fixo 
        #Condicao de Contorno
        self.sw_a = self.swc__ #gas  e espuma
        #self.sw_a = 1 # agua

    def create_circle_matrix(self,n):
    # n é o tamanho da matriz n x n
        matrix = np.zeros((n, n), dtype=int)
        center = 56  # O centro da matriz em 0.1mm
        radius = 56  # O maior raio possível em 0.1mm
        
        for i in range(n):
            for j in range(n):
                # Calcula a distância do centro
                distance = np.sqrt((i - center) ** 2 + (j - center) ** 2)
                # Define o valor como 1 se estiver dentro do círculo
                if distance <= radius:
                    matrix[i, j] = 1
                    
        return matrix

    # Exemplo de uso
   
    def settingEnviroment(self):

        self.tamMatrixEmXY = 113  # Tamanho da matriz tanto para X e Y em 0.1mm
        self.circle_matrix = self.create_circle_matrix(self.tamMatrixEmXY)
        self.h_x = 1 # distancia dos pontos da grid em y em 0.1mm
        self.h_y = 1 # distancia dos pontos da grid em y em 0.1mm
        self.h_t=0.1 #step de tempo em segundos
        self.timeT = 1 #tempo total para o deslocamento

        self.t = np.arange(0, self.timeT, self.h_t)

        self.tamX = self.tamMatrixEmXY #dimensão do sistema
        self.tamY = self.tamMatrixEmXY #dimensão do sistema
        self.steps = len(self.t) #numero de passos de tempo
        self.sol_tempo=[] # primeira fase (agua)
        self.sol_tempo2=[] # segunda fase (ar ou espuma)

        self.Sw = np.zeros([self.tamX,self.tamY])#x,y
        self.Sg = np.zeros([self.tamX,self.tamY])#x,y
        
        for i in range(self.tamX):
            for j in range(self.tamY):
                if self.circle_matrix[i,j] == 1:
                    self.Sw[i,j] = self.Sw_0
                else:
                    self.Sw[i,j] = -1

        self.sol_tempo.append(self.Sw)
        self.sol_tempo2.append(1-self.Sw)
        self.Sw_new=np.zeros([self.tamX,self.tamY])

    def calculoComDirecaoDoFluxo(self, vx,vy,fw,qx,qy,i,j):

        if vx>0 and vy>0:
            result = qx*(fw[i-1,j]  - fw[i,j] )  -qy*(fw[i,j-1]  - fw[i,j] )

        elif vx<0 and vy>0:
            result = qx*(fw[i+1,j]  - fw[i,j] )  -qy*(fw[i,j-1]  - fw[i,j] )

        elif vx>0 and vy<0:
            result = qx*(fw[i-1,j]  - fw[i,j] )  -qy*(fw[i,j+1]  - fw[i,j] )

        elif vx<0 and vy<0:
            result = qx*(fw[i+1,j]  - fw[i,j] )  -qy*(fw[i,j+1]  - fw[i,j] )

        
        return result

    def calculate(self):
        vx = 0.01294165212
        vy = 0.01294165212
        qx = (vx*self.h_t)/(self.h_x*self.phi__)
        qy = (vy*self.h_t)/(self.h_y*self.phi__)
        D = 2.299e-9 # agua à 25º é 2.299·10−9 m2·s−1
        fw = np.zeros([self.tamX,self.tamY])#x,y
        listOfFw = []
        self.Sw[0,56]=self.sw_a
        for i in range(self.tamX):
            for j in range(self.tamY):
                if self.circle_matrix[i,j] == 1:
                    if self.circle_matrix[i,j] == 1:
                            fw[i,j] = fw_(sw = self.Sw[i,j], sg = self.sg__, mrf = self.mrf__,krw0=self.krw0__,krg0=self.krg0__,lamb=self.lamb__,swc=self.swc__)
        
  
        listOfFw.append(np.copy(fw))
    
        for k in range(self.steps):
            for i in range(self.tamX):
                for j in range(self.tamY):
                    if self.circle_matrix[i,j] == 1:
                        qfw = 0
                        if i <112 and i>0 and j>0 and j<112:
                            if self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                # tem \/ and -> and <-
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                # tem /\ and -> and <-
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 1:
                                # tem \/ and <- and /\
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 1:
                                # tem /\ and -> and \/
                                ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                # tem \/ and ->
                                ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                # tem /\ and ->
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                # tem \/ and <-
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                # tem /\ and <-
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 0:
                                # tem <-
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x) + ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                # tem /\
                                ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 0:
                                # tem ->
                                ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 0 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                # tem \/
                                ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                            else:
                                #tem td
                                ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                        else:
                            if 0 == i:
                                if self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 1:
                                    # tem \/ and <- and /\
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                    # tem \/ and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                                elif self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                    # tem /\ and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i+1,j] == 1 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 0:
                                    # tem <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x) + ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                            elif 112 == i:
                                if self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 1:
                                    # tem \/ and -> and /\
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i,j-1] == 1 and self.circle_matrix[i,j+1] == 0:
                                    # tem \/ and ->
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)

                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 1:
                                    # tem /\ and ->
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i,j-1] == 0 and self.circle_matrix[i,j+1] == 0:
                                    # tem ->
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x) + ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                            elif 0 == j:
                                if self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 1  and self.circle_matrix[i,j+1] == 1:
                                    # tem /\ and -> and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1  and self.circle_matrix[i,j+1] == 1:
                                    # tem /\  and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0  and self.circle_matrix[i,j+1] == 1:
                                    # tem /\ and ->
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 0  and self.circle_matrix[i,j+1] == 1:
                                    # tem /\ 
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j+1] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                            elif 112 == j:
                                if self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 1  and self.circle_matrix[i,j] == 1:
                                    # tem \/ and -> and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 1  and self.circle_matrix[i,j] == 1:
                                    # tem \/  and <-
                                    ds = ((self.Sw[i+1,j] - 2*self.Sw[i,j] + self.Sw[i,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 1 and self.circle_matrix[i+1,j] == 0  and self.circle_matrix[i,j] == 1:
                                    # tem \/ and ->
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                                elif self.circle_matrix[i-1,j] == 0 and self.circle_matrix[i+1,j] == 0  and self.circle_matrix[i,j] == 1:
                                    # tem \/ 
                                    ds = ((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i-1,j])/self.h_x)+((self.Sw[i,j] - 2*self.Sw[i,j] + self.Sw[i,j-1])/self.h_y)
                                    qfw = self.calculoComDirecaoDoFluxo(vx,vy,fw,qx,qy,i,j)
                        # if(qfw==0):
                        #     print("deu 0")
                        # else:
                        #     print(self.Sw[i,j],qfw,ds)
                        self.Sw_new[i,j] = self.Sw[i,j] - qfw - D*ds


            self.Sw_new[0,56]=self.sw_a
            self.Sw = np.copy(self.Sw_new)
            self.sol_tempo.append(self.Sw)
            self.sol_tempo2.append(1-self.Sw)

            for i in range(0,self.tamX):
                for j in range(0,self.tamY):
                    if self.circle_matrix[i,j] == 1:
                        if self.circle_matrix[i,j] == 1:
                                fw[i,j] = fw_(self.Sw[i,j], self.sg__, self.mrf__,self.krw0__,self.krg0__,self.lamb__,self.swc__)
            listOfFw.append(np.copy(fw))        




    # saber se eles usam qual valor para fazer a conta do 

def fw_(sw,sg, mrf,krw0, krg0,lamb, swc):
    
    lambw = lambw_(sw =sw ,krw0 = krw0,lamb=lamb,swc=swc)
    lambt = lambw_(sw =sw ,krw0 = krw0,lamb=lamb,swc=swc)+lambg_(sg=sg,mrf=mrf,krg0=krg0,lamb=lamb,swc=swc,sw=sw)
    return lambw/lambt

def swe_(sw,swc,sgr):
    # if 1-swc-sgr == 0:
    #     print("deu 0")
    #     print(1,swc,sgr)
    #     asdfasdfasdf
    swe = 1
    if sw > swc and sw <= 1:
         swe = (sw-swc)/(1-swc-sgr)
    elif sw <= swc:
        Swe = 0
    else:
        Swe = 1
    return swe

def lambw_(sw,krw0,lamb,swc):
    muw = 1 # fixo
    krw = krw_(sw=sw,krw0=krw0,lamb=lamb,swc=swc)
    return (krw*sw)/muw

def lambg_(sg,mrf,krg0,lamb,swc,sw):
    mug = 0.0172 # fixo
    krg = krg_(krg0=krg0,lamb=lamb,swc=swc,sw=sw)
    return (krg*sg)/(mug*mrf)

def krw_(sw,krw0,lamb,swc):
    sgr = 0 # fixo 
    swe = swe_(sw,swc,sgr)
    return krw0*(swe**lamb)


def krg_(sw,krg0,lamb,swc):
    sgr = 0 # fixo 
    swe = swe_(sw,swc,sgr)
    if swe != 0: #por enquanto
        return krg0*(1-swe**(3-(2/lamb)))
    else:
        return krg0

def phi_(volW,volT):
    return volW/volT

