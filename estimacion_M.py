# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:16:23 2019

@author: carlo
"""
from gurobipy import *
import numpy as np
import neighborhood as neigh
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import auxiliar_functions as af


# def estima_BigM(data):
#
#     m = len(data)
#     BigM = 0
#
#     for i in range(m):
#         for j in range(m):
#             if i != j:
#                 comp1 = data[i]
#                 comp2 = data[j]
#
#                 if type(comp1) is neigh.Poligono or neigh.Poligonal:
#                     if type(comp2) is neigh.Poligono:
#                         maximo = max([np.linalg.norm(v-w) for v in comp1.V
#                                                           for w in comp2.V])
#
#                     if type(comp2) is neigh.Elipse:
#                         maximo = comp2.radii + max([np.linalg.norm(v-comp2.centfro) for v in comp1.V])
#
#                 if type(comp1) is neigh.Poligonal:
#                     if type(comp2) is neigh.Poligono or neigh.Poligonal:
#                         maximo = max([np.linalg.norm(w-v) for w in comp2.V
#                                                           for v in comp1.V])
#
#                     if type(comp2) is neigh.Elipse:
#                         caso1 = comp2.radii + np.linalg.norm(comp1.P-comp2.centro)
#                         caso2 = comp2.radii + np.linalg.norm(comp1.Q-comp2.centro)
#                         maximo = max(caso1, caso2)
#
#                 if type(comp1) is neigh.Elipse:
#                     if type(comp2) is neigh.Poligono or neigh.Poligonal:
#                         maximo = comp1.radii + max([np.linalg.norm(comp1.centro-w) for w in comp2.V])
#
#                     if type(comp2) is neigh.Elipse:
#                         maximo = comp1.radii + np.linalg.norm(comp1.centro-comp2.centro) + comp2.radii
#
#                 if maximo >= BigM:
#                     BigM = maximo
#
#     return BigM


def estima_BigM_local(comp1, comp2):
        maximo = 0
        if type(comp1) is neigh.Poligono or type(comp1) is neigh.Poligonal:
            if type(comp2) is neigh.Poligono or type(comp2) is neigh.Poligonal:
                maximo = max([np.linalg.norm(v-w) for v in comp1.V
                                                  for w in comp2.V])

            if type(comp2) is neigh.Elipse:
                maximo = comp2.radii + max([np.linalg.norm(v-comp2.centro) for v in comp1.V])

        if type(comp1) is neigh.Elipse:
            if type(comp2) is neigh.Poligono or type(comp2) is neigh.Poligonal:
                maximo = comp1.radii + max([np.linalg.norm(comp1.centro-w) for w in comp2.V])

            if type(comp2) is neigh.Elipse:
                maximo = comp1.radii + np.linalg.norm(comp1.centro-comp2.centro) + comp2.radii

        return maximo

def estima_SmallM_local(comp1, comp2):
        if type(comp1) is neigh.Poligono or type(comp1) is neigh.Poligonal:
            if type(comp2) is neigh.Poligono or type(comp2) is neigh.Poligonal:
                minimo = min([np.linalg.norm(v-w) for v in comp1.V
                                                  for w in comp2.V])

            if type(comp2) is neigh.Elipse:
                minimo = - comp2.radii + min([np.linalg.norm(v-comp2.centro) for v in comp1.V])

        if type(comp1) is neigh.Elipse:
            if type(comp2) is neigh.Poligono or type(comp2) is neigh.Poligonal:
                minimo = -comp1.radii + min([np.linalg.norm(comp1.centro-w) for w in comp2.V])

            if type(comp2) is neigh.Elipse:
                minimo = -comp1.radii + np.linalg.norm(comp1.centro-comp2.centro) - comp2.radii

        return minimo

def estima_max_inside(comp):
        maximo = 0
        if type(comp) is neigh.Poligono:
            maximo = max([np.linalg.norm(v-w) for v in comp.V for w in comp.V])

        if type(comp) is neigh.Elipse:
            maximo = 2*comp.radii

        if type(comp) is neigh.Poligonal:
            maximo = comp.alpha * comp.longitud

        return maximo

def estima_M_alpha1(entorno, punto1, punto2):
    if type(entorno) is neigh.Circle:

        theta = np.linspace(0, 2*np.pi, 20)

        centro = entorno.center
        radii = entorno.radii

        x = centro[0] + radii*np.cos(theta)
        y = centro[1] + radii*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant([x[i], y[i]], punto1, punto2) for i in range(20)]



    if type(entorno) is neigh.Poligonal:

        theta = np.linspace(0, 1, 20)

        extreme_point = np.array(entorno.V[0])

        vector = np.array(entorno.V[1]) - np.array(entorno.V[0])

        determinantes = [af.determinant(extreme_point + theta[i]*vector, punto1, punto2) for i in range(20)]

    m = min(determinantes)/2
    M = 2*max(determinantes)

    return m, M
def estima_M_alpha2(punto1, entorno, punto2):
    if type(entorno) is neigh.Circle:

        theta = np.linspace(0, 2*np.pi, 20)

        centro = entorno.center
        radii = entorno.radii

        x = centro[0] + radii*np.cos(theta)
        y = centro[1] + radii*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, [x[i], y[i]], punto2) for i in range(20)]

    if type(entorno) is neigh.Poligonal:
        theta = np.linspace(0, 1, 20)

        extreme_point = np.array(entorno.V[0])

        vector = np.array(entorno.V[1]) - np.array(entorno.V[0])

        determinantes = [af.determinant(punto1, extreme_point + theta[i]*vector, punto2) for i in range(20)]

    m = min(determinantes)/2
    M = 2*max(determinantes)

    return m, M

def estima_M_alpha3(punto1, punto2, entorno):
    if type(entorno) is neigh.Circle:

        theta = np.linspace(0, 2*np.pi, 20)

        centro = entorno.center
        radii = entorno.radii

        x = centro[0] + radii*np.cos(theta)
        y = centro[1] + radii*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, punto2, [x[i], y[i]]) for i in range(20)]

    if type(entorno) is neigh.Poligonal:
        theta = np.linspace(0, 1, 20)

        extreme_point = np.array(entorno.V[0])

        vector = np.array(entorno.V[1]) - np.array(entorno.V[0])

        determinantes = [af.determinant(punto1, punto2, extreme_point + theta[i]*vector) for i in range(20)]

    m = min(determinantes)/2
    M = 2*max(determinantes)

    return m, M

def estima_M_alpha4(punto1, entorno1, entorno2):
    if type(entorno1) is neigh.Circle and type(entorno2) is neigh.Circle:

        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = entorno1.center
        radii1 = entorno1.radii

        x1 = centro1[0] + radii1*np.cos(theta)
        y1 = centro1[1] + radii1*np.sin(theta)

        centro2 = entorno2.center
        radii2 = entorno2.radii

        x2 = centro2[0] + radii2*np.cos(theta)
        y2 = centro2[1] + radii2*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, [x1[i], y1[i]], [x2[j], y2[j]]) for i in range(20) for j in range(20)]

    if type(entorno1) is neigh.Poligonal and type(entorno2) is neigh.Poligonal:
        theta = np.linspace(0, 1, 20)

        extreme_point1 = np.array(entorno1.V[0])

        vector1 = np.array(entorno1.V[1]) - np.array(entorno1.V[0])

        extreme_point2 = np.array(entorno2.V[0])

        vector2 = np.array(entorno2.V[1]) - np.array(entorno2.V[0])

        determinantes = [af.determinant(punto1, extreme_point1 + theta[i]*vector1, extreme_point2 + theta[j]*vector2) for i in range(20) for j in range(20)]

    m = min(determinantes)/2
    M = 2*max(determinantes)

    return m, M

# def estima_L(neighborhood, punto):
#     if type(neighborhood) is neigh.Circle:
#         centro = neighborhood.center
#         dis = np.linalg.norm(np.array(centro) - np.array(punto)) - neighborhood.radii
#
#         # dist = 0
#
#     if type(neighborhood) is neigh.Poligonal:
#         dis = af.dist_point_segment(punto, neighborhood)
#         # dr = np.array(neighborhood.V[1]) - np.array(neighborhood.V[0])
#         #
#         # A, B = -dr[1], dr[0]
#         #
#         # C = -A*neighborhood.V[0][0] - B*neighborhood.V[0][1]
#         #
#         # dis = abs(A*punto[0] + B*punto[1] + C)/np.sqrt(A**2 + B**2)
#
#     return dis
#
#
# def estima_U(neighborhood, punto):
#     if type(punto) is neigh.Circle:
#         centro1 = neighborhood.center
#         centro2 = punto.center
#
#         dis = np.linalg.norm(np.array(centro1) - np.array(centro2)) + neighborhood.radii + punto.radii
#
#     else:
#         if type(neighborhood) is neigh.Circle:
#             dis = np.linalg.norm(np.array(neighborhood.center) - np.array(punto)) + neighborhood.radii
#
#             # dist = 10000
#
#         if type(neighborhood) is neigh.Poligonal:
#             dis = max([np.linalg.norm(np.array(neighborhood.V[i]) - np.array(punto)) for i in range(2)])
#
#     return dis


def estima_M_complete(ent1, ent2):

    if type(ent1) is neigh.Elipse and type(ent2) is neigh.Elipse:
        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        centro2 = ent2.center
        radii21 = ent2.width
        radii22 = ent2.height

        x2 = centro2[0] + radii21*np.cos(theta)
        y2 = centro2[1] + radii22*np.sin(theta)

        # print([x[0], y[0]])

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(20) for j in range(20)]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if type(ent1) is neigh.Elipse and type(ent2) is neigh.Circle:
        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        centro2 = ent2.center
        radii2 = ent2.radii

        x2 = centro2[0] + radii2*np.cos(theta)
        y2 = centro2[1] + radii2*np.sin(theta)

        # print([x[0], y[0]])

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(20) for j in range(20)]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if type(ent1) is neigh.Elipse and type(ent2) is neigh.Punto:
        theta = np.linspace(0, 2*np.pi, 20)

        centro = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x = centro[0] + radii11*np.cos(theta)
        y = centro[1] + radii12*np.sin(theta)

        # print([x[0], y[0]])

        distancias = [np.linalg.norm(np.array([x[i], y[i]]) - np.array(ent2.V)) for i in range(20)]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if type(ent1) is neigh.Elipse and (type(ent2) is neigh.Poligono or type(ent2) is neigh.Poligonal):

        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array(v)) for i in range(20) for v in ent2.V]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if type(ent1) is neigh.Circle and type(ent2) is neigh.Elipse:
        return estima_M_complete(ent2, ent1)

    if type(ent1) is neigh.Circle and type(ent2) is neigh.Circle:
        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = ent1.center
        radii1 = ent1.radii

        x1 = centro1[0] + radii1*np.cos(theta)
        y1 = centro1[1] + radii1*np.sin(theta)

        centro2 = ent2.center
        radii2 = ent2.radii

        x2 = centro2[0] + radii2*np.cos(theta)
        y2 = centro2[1] + radii2*np.sin(theta)

        # print([x[0], y[0]])

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(20) for j in range(20)]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if type(ent1) is neigh.Circle and type(ent2) is neigh.Punto:
        theta = np.linspace(0, 2*np.pi, 20)

        centro = ent1.center
        radii1 = ent1.radii

        x = centro[0] + radii1*np.cos(theta)
        y = centro[1] + radii1*np.sin(theta)

        # print([x[0], y[0]])

        distancias = [np.linalg.norm(np.array(ent2.V) - np.array([x[i], y[i]])) for i in range(20)]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U


    if type(ent1) is neigh.Circle and (type(ent2) is neigh.Poligono or type(ent2) is neigh.Poligonal):
        theta = np.linspace(0, 2*np.pi, 20)

        centro1 = ent1.center
        radii1 = ent1.radii

        x1 = centro1[0] + radii1*np.cos(theta)
        y1 = centro1[1] + radii1*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array(v)) for i in range(20) for v in ent2.V]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U


    if type(ent1) is neigh.Punto and type(ent2) is neigh.Elipse:
        return estima_M_complete(ent2, ent1)

    if type(ent1) is neigh.Punto and type(ent2) is neigh.Circle:
        return estima_M_complete(ent2, ent1)

    if type(ent1) is neigh.Punto and type(ent2) is neigh.Punto:
        L = np.linalg.norm(np.array(ent1.V) - np.array(ent2.V))/2
        U = 2*L
        return L, U

    if type(ent1) is neigh.Punto and (type(ent2) is neigh.Poligono or type(ent2) is neigh.Poligonal):
        distancias = [np.linalg.norm(np.array(ent1.V) - np.array(v)) for v in ent2.V]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U

    if (type(ent1) is neigh.Poligono or type(ent1) is neigh.Poligonal) and type(ent2) is neigh.Elipse:
        return estima_M_complete(ent2, ent1)

    if (type(ent1) is neigh.Poligono or type(ent1) is neigh.Poligonal) and type(ent2) is neigh.Circle:
        return estima_M_complete(ent2, ent1)

    if (type(ent1) is neigh.Poligono or type(ent1) is neigh.Poligonal) and type(ent2) is neigh.Punto:
        return estima_M_complete(ent2, ent1)

    if (type(ent1) is neigh.Poligono or type(ent1) is neigh.Poligonal) and (type(ent2) is neigh.Poligono or type(ent2) is neigh.Poligonal):
        distancias = [np.linalg.norm(np.array(v) - np.array(w)) for v in ent1.V for w in ent2.V]

        L = min(distancias)/2
        U = 2*max(distancias)

        return L, U


