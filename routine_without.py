import gurobipy as gp
import pdb
# from HTSPS_without_prepro import HTSPS_without_prepro
from neighborhood import Circle
from tspn_b import tspn_b

import numpy as np
import pandas as pd

init = False

# if init:
#     dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
# else:
#     dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])

dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])

for nP in [5, 10, 20, 30, 50, 80, 90, 100]:
    for instance in range(10):
        
        print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')
        
        segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        barriers = []
        for lista in segments:
            barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
            
        bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        neighborhoods = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

        resultados = tspn_b(barriers, neighborhoods, prepro=False, log=False, dominant = False, picture=False, time_limit=3600, init=False)
        
        serie = pd.Series([instance] + resultados, index = dataframe.columns)
        
        dataframe = dataframe.append(serie, ignore_index=True)
        dataframe.to_csv('./resultados/resultados_without_prepro_bolas_visitar_updated.csv')

        
