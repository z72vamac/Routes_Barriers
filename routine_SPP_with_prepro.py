import gurobipy as gp
import pdb
from sppn_b import HSPPS
from neighborhood import Circle

import numpy as np
import pandas as pd

init = False

if init:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
else:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'NodeCount', 'ObjVal'])

for nP in [5, 10, 20, 30, 50, 80]:
    for instance in range(10):
        
        print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')
        
        segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        barriers = []
        for lista in segments:
            barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
            
        bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        N = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

        resultados = HSPPS(barriers, N, timeLimit = 3600)
        
        serie = pd.Series([instance] + resultados, index = dataframe.columns)
        
        dataframe = dataframe.append(serie, ignore_index=True)
        dataframe.to_csv('./resultados/SPP_with_prepro.csv')

        
