import gurobipy as gp
import pdb
# from HTSPS_with_prepro import HTSPS_with_prepro
from neighborhood import Circle, Poligonal
from tspn_b import tspn_b

import numpy as np
import pandas as pd

init = False

if init:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
else:
    dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal'])

dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])

for nP in [80, 90, 100]:
#for nP in [5, 10, 20, 30, 50, 80, 90, 100]:

    for instance in range(10):
        
        print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')
        
        segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        barriers = []
        for lista in segments:
            barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])
            
        segmentos_visitar = np.genfromtxt('./instancias/segmentos_visitar' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

        neighborhoods = [Poligonal(V = [np.array([lista[0], lista[1]]), np.array([lista[2], lista[3]])]) for lista in segmentos_visitar]

        resultados = tspn_b(barriers, neighborhoods, prepro=False, log=False, dominant = False, picture=False, time_limit=3600, init=False)
        
        serie = pd.Series([instance] + resultados, index = dataframe.columns)
        
        dataframe = dataframe.append(serie, ignore_index=True)
        dataframe.to_csv('./resultados/resultados_without_prepro_segmentos_visitar_updated-2.csv')

        
