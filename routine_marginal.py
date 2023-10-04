# Routine that evaluates how each strengthening step affects to the solver behaviour.

"""
Computing the routine for Smith instances
"""

import gurobipy as gp
import pdb
# from HTSPS_with_prepro import HTSPS_with_prepro
from neighborhood import Circle
from tspn_b import tspn_b

import numpy as np
import pandas as pd

init = False

dataframe = pd.DataFrame(columns=['Instance', 'Radii', 'n_N', 'n_B', 'A4', 'edge_removing', 'variable_fixing', 'bigM_estimation', 'Gap', 'Runtime', 'Time_Prepro', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])


start = False

num_rows = 0
if start:
    dataframe = pd.read_csv('./resultados/results_smith_marginal.csv').iloc[:, 1:]
    num_rows = dataframe.shape[0]-1

counter = 1

for nP in [10]:
    for instance in range(10):
        for radii in [1]:
            if counter > num_rows:
                print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')

                segments = np.genfromtxt('./instancias_smith/barreras'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                print('A4 = ' + str(False))
                print('Removing edges: ' + str(True))
                print('Fixing variables: ' + str(False))
                print('Estimating BigM: ' + str(False))
                
                barriers = []
                for lista in segments:
                    barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                nB = len(barriers)

                bolas = np.genfromtxt('./instancias_smith/bolas'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                neighborhoods = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

                resultados = tspn_b(barriers, neighborhoods, A4=False, edge_removing=True, variable_fixing=False, bigM_estimation=False, log=False, dominant = False, picture=False, time_limit=3600, init=False)

                serie = pd.Series([instance, radii] + resultados, index = dataframe.columns)

                dataframe = dataframe.append(serie, ignore_index=True)
                dataframe.to_csv('./resultados/results_smith_marginal.csv')

            counter +=1

            if counter > num_rows:
                print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')

                segments = np.genfromtxt('./instancias_smith/barreras'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                print('A4 = ' + str(False))
                print('Removing edges: ' + str(False))
                print('Fixing variables: ' + str(True))
                print('Estimating BigM: ' + str(False))
                
                barriers = []
                for lista in segments:
                    barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                nB = len(barriers)

                bolas = np.genfromtxt('./instancias_smith/bolas'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                neighborhoods = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

                resultados = tspn_b(barriers, neighborhoods, A4=False, edge_removing=False, variable_fixing=True, bigM_estimation=False, log=False, dominant = False, picture=False, time_limit=3600, init=False)

                serie = pd.Series([instance, radii] + resultados, index = dataframe.columns)

                dataframe = dataframe.append(serie, ignore_index=True)
                dataframe.to_csv('./resultados/results_smith_marginal.csv')

            counter += 1

            if counter > num_rows:
                print('\n\nResolviendo la instancia ' + str(instance) + ' con un numero ' + str(nP) + ' de neighborhoods.\n\n')

                segments = np.genfromtxt('./instancias_smith/barreras'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                print('A4 = ' + str(False))
                print('Removing edges: ' + str(False))
                print('Fixing variables: ' + str(False))
                print('Estimating BigM: ' + str(True))
                
                barriers = []
                for lista in segments:
                    barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                nB = len(barriers)

                bolas = np.genfromtxt('./instancias_smith/bolas'+ str(nP) + '-' + str(instance) + '-' + str(radii) + '.csv', delimiter = ',')

                neighborhoods = [Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

                resultados = tspn_b(barriers, neighborhoods, A4=False, edge_removing=False, variable_fixing=False, bigM_estimation=True, log=False, dominant = False, picture=False, time_limit=3600, init=False)

                serie = pd.Series([instance, radii] + resultados, index = dataframe.columns)

                dataframe = dataframe.append(serie, ignore_index=True)
                dataframe.to_csv('./resultados/results_smith_marginal.csv')

            counter += 1

