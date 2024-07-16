import numpy as np
import math
import matplotlib.pyplot as plt
from codigo import PorositiyProblem

if __name__ == '__main__':
    pp = PorositiyProblem()
    pp.settingEnviroment()

    # plt.matshow(pp.Sw)
    # plt.show()
    solcmep = np.arange(0,1.001,0.001)
    pp.calculate()
    for i in range(0,len(pp.sol_tempo)):
        plt.matshow(pp.sol_tempo2[i] )
        plt.show()


