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



# segments = np.genfromtxt('./instancias/segmentos30-2.csv', delimiter = ',')

# barriers = []
# for lista in segments[25:28]:
#     barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

# bolas = np.genfromtxt('./instancias/bolas30-2.csv', delimiter = ',')[2:4]
# N = [neigh.Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]
# segmentos_visitar = np.genfromtxt('./instancias/segmentos_visitar30-2.csv', delimiter = ',')
# N = [neigh.Poligonal(V = [np.array([lista[0], lista[1]]), np.array([lista[2], lista[3]])]) for lista in segmentos_visitar] # 105.164


# Example of Figures 2 and 3
figure = 2 # or 3

segments = np.genfromtxt('./instance_example/barriers{0}.csv'.format(figure), delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

circles = np.genfromtxt('./instance_example/circles.csv'.format(figure), delimiter = ',')

neighbourhoods = [neigh.Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]

# Random instances
instance = 5
neighbourhood_size = 10

segments = np.genfromtxt('./instances_random/barriers/barriers{0}-{1}.csv'.format(neighbourhood_size, instance), delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

circles = np.genfromtxt('./instances_random/circles/circles{0}-{1}.csv'.format(neighbourhood_size, instance), delimiter = ',')

neighbourhoods = [neigh.Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]

resultados = tspn_b(barriers, neighbourhoods, A4 = False, dominant=False, log=False, picture=True, time_limit=600, init = False)
