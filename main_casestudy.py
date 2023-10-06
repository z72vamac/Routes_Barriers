import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
# from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
import neighborhood as neigh
import copy
import estimacion_M as eM
import auxiliar_functions as af
import networkx as nx
from HTSPS_with_prepro2 import HTSPS_with_prepro2
from HTSPS_without_prepro import HTSPS_without_prepro
# from HTSPS_with_prepro import HTSPS_with_prepro
# from HTSPS_with_prepro3 import HTSPS_with_prepro3
# from HTSPS_with_prepro4 import HTSPS_with_prepro4
# from HTSPS_new_ven import tspn_b
from tspn_b import tspn_b
from tspn_b2 import tspn_b2
from heuristic import heuristic
from heuristic2 import heuristic2

# from sppn_b import sppn_b
from HTSPS_new_ven import HTSPS_ven

# from HTSPS_without_prepro import HTSPS_without_prepro


blocks = [2]
barriers = []

for i in blocks:

    segments = np.genfromtxt('./case_study/barriers{0}.csv'.format(i), delimiter = ',')

    for lista in segments:
        barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

print(len(barriers)) 
bolas = [[59.5, 57.5, 7.5], [(57+65.5)/2, 200-(92+100.5)/2, 4.25], [(77+81)/2, 200-82, 2], [(98+103)/2, 200-(85+90)/2, 2.5], [(117+124)/2, 200-(100+107)/2, 3.5]]

N = [neigh.Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

resultados = tspn_b(barriers, N, A4 = False, dominant=False, prepro=True, log=False, picture=True, time_limit=3*3600, init = False)


# print(resultados)

# HTSPS_with_prepro(barriers, N)
