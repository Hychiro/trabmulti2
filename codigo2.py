import numpy as np
import math
import matplotlib.pyplot as plt

class PorositiyProblem():
    def __init__(self,krw0 = 0.7737, krg =  0.2007, lamb = 0.52, swc = 0.1, mrf =1.0):

        self._krw0 = krw0 # ajustavel
        self._krg0 = krg # ajustavel
        self._lamb = lamb # ajustavel
        self._swc = swc # ajustavel
        self._mrf = mrf # valor referente ao sufactante, ajustavel se n for agua

        self._sw = 1.0 # vai variar na conta, mas é o valor inicial da saturação
        self._sg = 0.0 # constante, sg é o resultado do sw-swc se me lembro bem (saturação de liquido - saturação conata = saturação do gas) já q sgr pode ser considerado 0
        

        self._volW = 110.893 # volume de agua dentro do experimento (resultado da massaFinal - massaInicial do experimento)
        self._volt = 319.785 # volume total do ambiente
        self._phi = self.phi_() # valor da porosidade
        
        self.tamMatrixEmXY = 113  # Tamanho da matriz tanto para X e Y em 1mm
        
        #Condicao inicial
        self._Sw_0 = 1.0 #fixo 
        #Condicao de Contorno
        # self._sw_a = self._swc #gas  e espuma
        self._sw_a = 1 # agua

    def create_circle_matrix(self,n):
    # n é o tamanho da matriz n x n
        matrix = np.zeros((n, n), dtype=int)
        center = 56  # O centro da matriz em 1mm
        radius = 56  # O maior raio possível em 1mm
        
        for i in range(n):
            for j in range(n):
                # Calcula a distância do centro
                distance = np.sqrt((i - center) ** 2 + (j - center) ** 2)
                # Define o valor como 1 se estiver dentro do círculo
                if distance <= radius:
                    matrix[i, j] = 1
                    
        return matrix
    
    # Exemplo de uso
   
    def settingEnviroment(self,timeT = 20):

        self.circle_matrix = self.create_circle_matrix(self.tamMatrixEmXY)

        self.h_x = 1 # distancia dos pontos da grid em y em 1mm
        self.h_y = 1 # distancia dos pontos da grid em y em 1mm
        self.h_t = 0.1 #step de tempo em segundos
        self._timeT = timeT #tempo total para o deslocamento

        self.t = np.arange(0, self._timeT, self.h_t)

        self.tamX = self.tamMatrixEmXY #dimensão do sistema
        self.tamY = self.tamMatrixEmXY #dimensão do sistema

        
        self.sol_tempo=[] # primeira fase (agua)
        self.sol_tempo2=[] # segunda fase (ar ou espuma)

        self.Sw = np.zeros([len(self.t),self.tamX,self.tamY])
        # self.Sw = np.zeros([self.tamX,self.tamY])#x,y
        # self.Sg = np.zeros([self.tamX,self.tamY])#x,y
        
        for i in range(self.tamX):
            for j in range(self.tamY):
                self.Sw[0,i,j] = 1

        # self.sol_tempo.append(self.Sw)
        # self.sol_tempo2.append(1-self.Sw)
        self.Sw_new=np.zeros([self.tamX,self.tamY])


        # saber se eles usam qual valor para fazer a conta do 


    def calculate(self):
        vx = 0.01294165212 #velociade x em m/s ?
        vy = 0.01294165212 #velocidade y em m/s ?
        # qx = (vx*self.h_t)/(self.h_x*self._phi)
        # qy = (vy*self.h_t)/(self.h_y*self._phi)
        D = 2.299e-3 # constante de difussão da agua à 25º é 2.299·10−9 m2·s−1
        fw = np.zeros([len(self.t),self.tamX,self.tamY])
        for i in range(self.tamX):
            for j in range(self.tamY):
                if self.circle_matrix[i,j] == 1:
                    fw[0,i,j] = self.fw_(sw = self.Sw[0,i,j], sg = self._sg, mrf = self._mrf,krw0=self._krw0,krg0=self._krg0,lamb=self._lamb,swc=self._swc)

        for k in range(1,len(self.t)):
            for i in range(1,self.tamX-1):
                for j in range(1,self.tamY-1):
                    if self.circle_matrix[i,j] == 1:
                        ds = ((D *self.h_t)/self._phi)((self.Sw[k-1,i+1,j] - 2*self.Sw[k-1,i,j] + self.Sw[k-1,i-1,j])/(self.h_x**2))+((self.Sw[k-1,i,j+1] - 2*self.Sw[k-1,i,j] + self.Sw[k-1,i,j-1])/(self.h_y**2))
                        qfw = -(((vx * self.h_t)/self._phi) * (fw[k-1,i+1,j]-fw[k-1,i-1,j])/(2*self.h_x)) - (((vy * self.h_t)/self._phi)*(fw[k-1,i,j+1]-fw[k-1,i,j-1])/(2*self.h_y))
                        self.Sw[k,i,j] = self.Sw[k-1,i,j] + qfw + ds

            for i2 in range(self.tamX):
                for j2 in range(self.tamY):
                    if self.circle_matrix[i2,j2] == 1:
                        fw[k,i2,j2] = self.fw_(sw = self.Sw[k,i2,j2], sg = self._sg, mrf = self._mrf,krw0=self._krw0,krg0=self._krg0,lamb=self._lamb,swc=self._swc)
    
    def fw_(self,sw,sg, mrf,krw0, krg0,lamb, swc):
        
        lambw = self.lambw_(sw =sw ,krw0 = krw0,lamb=lamb,swc=swc)
        lambt = self.lambw_(sw =sw ,krw0 = krw0,lamb=lamb,swc=swc)+self.lambg_(sg=sg,mrf=mrf,krg0=krg0,lamb=lamb,swc=swc,sw=sw)
        return lambw/lambt

    def swe_(self,sw,swc,sgr):
        # if 1-swc-sgr == 0:
        #     print("deu 0")
        #     print(1,swc,sgr)
        #     asdfasdfasdf
        swe = 1
        if sw > swc and sw <= 1:
            swe = (sw-swc)/(1-swc-sgr)
        elif sw <= swc:
            swe = 0
        else:
            swe = 1
        return swe

    def lambw_(self,sw,krw0,lamb,swc):
        muw = 1 # fixo
        krw = self.krw_(sw=sw,krw0=krw0,lamb=lamb,swc=swc)
        return (krw*sw)/muw

    def lambg_(self,sg,mrf,krg0,lamb,swc,sw):
        mug = 0.0172 # fixo
        krg = self.krg_(krg0=krg0,lamb=lamb,swc=swc,sw=sw)
        return (krg*sg)/(mug*mrf)

    def krw_(self,sw,krw0,lamb,swc):
        sgr = 0 # fixo 
        swe = self.swe_(sw,swc,sgr)
        return krw0*(swe**lamb)


    def krg_(self,sw,krg0,lamb,swc):
        sgr = 0 # fixo 
        swe = self.swe_(sw,swc,sgr)
        if swe != 0: #por enquanto
            return krg0*(1-swe**(3-(2/lamb)))
        else:
            return krg0

    def phi_(self):
        return self._volW/self._volt


