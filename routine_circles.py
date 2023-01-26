import gurobipy as gp
import pdb
# from HTSPS_with_prepro import HTSPS_with_prepro
from neighborhood import Circle
from tspn_b import tspn_b

import numpy as np
import pandas as pd

init = False

# if init:
#     dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Strength', 'A4', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])
# else:
#     dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Strength', 'A4', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal'])

dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'Strength', 'A4', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])

A4s = [False, True]
prepros = [False, True]

start = True

num_rows = 0
if start:
    dataframe = pd.read_csv('./resultados/results_circles_70.csv').iloc[:, 1:]
    num_rows = dataframe.shape[0] - 1

counter = 1

for nP in [70, 75, 80]: #[5, 10, 20, 30, 50, 80, 100]:
    for a4 in A4s:
        for prepro in prepros:
            for instance in range(5):
                if counter > num_rows:
                    print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')

                    segments = np.genfromtxt('./instancias/segmentos'+ str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

                    print('A4 = ' + str(a4))
                    print('Strengthening: ' + str(prepro))

                    barriers = []
                    for lista in segments:
                        barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                    nB = len(barriers)

                    if not(a4):
                        sublist = np.random.choice(nB, int(np.floor(0.5 * nB)))
                        barriers = [barriers[b] for b in sublist]

                    bolas = np.genfromtxt('./instancias/bolas' + str(nP) + '-' + str(instance) + '.csv', delimiter = ',')

                    neighborhoods = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

                    resultados = tspn_b(barriers, neighborhoods, prepro=prepro, A4 = a4, log=False, dominant = False, picture=False, time_limit=3600, init=False)

                    serie = pd.Series([instance] + resultados, index = dataframe.columns)

                    dataframe = dataframe.append(serie, ignore_index=True)
                    dataframe.to_csv('./resultados/results_circles_70.csv')

                counter += 1

        
