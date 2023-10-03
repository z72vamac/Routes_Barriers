# TSPN-B
# Resolviendo el TSPN-B obligando a que los puntos sean el centro.

import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh


def heuristic2(barriers, neighborhoods, prepro=True, log=False, dominant=False, picture=False, time_limit=7200, init=False):

    first_time = time.time()

    vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
    vertices_neighborhood = vertices_neighborhood[::-1]

    edges_neighborhood = []

    for (a, b) in vertices_neighborhood:
        for c in range(len(barriers)):
            for d in range(2):
                # if prepro:
                    # Point of the barrier to check if is visible by the neighborhood
                point = barriers[c][d]

                # Neighborhood to check if it is visible by the point
                neighborhood = neighborhoods[abs(a) - 1]

                if type(neighborhood) is neigh.Circle:
                    center = neighborhood.center

                if type(neighborhood) is neigh.Poligonal:
                    center = np.mean(neighborhood.V, axis = 1)

                barrier = [[point[0], point[1]], [center[0], center[1]]]

                intersect = False
                for barrieri in barriers:
                    if af.intersect(barrieri, barrier):
                        intersect = True
                        break

                if not(intersect):
                    edges_neighborhood.append((a, b, c, d))
                    edges_neighborhood.append((c, d, a, b))

    vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))

    edges_barrier = []
    for v, i in vertices_barrier:
        for w, j in vertices_barrier:
            if v < w:
                # if prepro:
                barrier = [barriers[v][i], barriers[w][j]]
                if np.linalg.norm(np.array(barriers[v][i]) - np.array(barriers[w][j])) >= 0.5:
                    intersect = False
                    for barrieri in barriers:
                        if af.intersect(barrieri, barrier):
                            intersect = True
                            break

                    if not (intersect):
                        edges_barrier.append((v, i, w, j))
                        edges_barrier.append((w, j, v, i))

                # else:
                #     edges_barrier.append((v, i, w, j))

    indices_barriers = [(v, 0, v, 1) for v in range(len(barriers))]

    vertices_total = vertices_neighborhood + vertices_barrier
    edges_total = edges_neighborhood + edges_barrier

    if log:
        print("vertices_neighborhood = " + str(vertices_neighborhood))
        print("vertices_barrier = " + str(vertices_barrier))

        print("edges_neighborhood = " + str(edges_neighborhood))
        print("edges_barrier = " + str(edges_barrier))

    y_index = edges_total

    if log:
        print("y_index = " + str(y_index))

    # P_S and P_T: indices of the points in the neighborhoods
    # p_index = []
    # for a, b, c, d in edges_total:
    #     for e, f in vertices_neighborhood:
    #         p_index.append((a, b, c, d, e))
    p_index = edges_neighborhood
    # for index in vertices_neighborhood:
    #     for dim in range(2):
    #         p_index.append((index[0], index[1], dim))

    if log:
        print("p_index = " + str(p_index))

    dist = {}

    for a, b, c, d in edges_total:
        # dist constraints
        if (a, b, c, d) in edges_barrier:
            dist[a, b, c, d] = np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d]))

        if (a, b, c, d) in edges_neighborhood:
            if a < 0:
                neighborhood = neighborhoods[abs(a) - 1]

                if type(neighborhood) is neigh.Circle:
                    center = neighborhood.center

                if type(neighborhood) is neigh.Poligonal:
                    center = np.mean(neighborhood.V, axis=1)

                dist[a, b, c, d] = np.linalg.norm(center - np.array(barriers[c][d]))

            if c < 0:
                neighborhood = neighborhoods[abs(c) - 1]

                if type(neighborhood) is neigh.Circle:
                    center = neighborhood.center

                if type(neighborhood) is neigh.Poligonal:
                    center = np.mean(neighborhood.V, axis=1)

                dist[a, b, c, d] = np.linalg.norm(center - np.array(barriers[a][b]))

    # if log:
    #     print("dist = " + str(dist_index))

    # f variables:
    g_index = y_index

    if log:
        print("g_index = " + str(g_index))

    model = gp.Model('HTSPS_Model')

    p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
    y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
    # dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')

    # z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    g_var = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')

    model.update()



        # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)

    # model.addConstrs(z[v, v] == 0 for v in vertices_total

    # Restriccion 1
    for v_n, j in vertices_neighborhood:
        model.addConstr(
            gp.quicksum(y[v, i, v_n, j] for v, i in vertices_total if (v, i, v_n, j) in edges_neighborhood) >= 1)

    # Restriccion 2
    for v, i in vertices_total:
        model.addConstr(gp.quicksum(y[v, i, v_p, j] for v_p, j in vertices_total if (v, i, v_p, j) in edges_total)
                        == gp.quicksum(
            y[v_p, j, v, i] for v_p, j in vertices_total if (v_p, j, v, i) in edges_total))

    # Restriccion 3
    for v_n, i in vertices_neighborhood:
        if v_n <= -2:
            model.addConstr(gp.quicksum(g_var[v_n, i, v, j] for v, j in vertices_total if
                                        (v_n, i, v, j) in edges_neighborhood) - gp.quicksum(
                g_var[v, j, v_n, i] for v, j in vertices_total if
                (v, j, v_n, i) in edges_neighborhood) == 1)

    # Restriccion 4
    for v_b, i in vertices_barrier:
        model.addConstr(gp.quicksum(g_var[(w, j, v_b, i)] for w, j in vertices_total if
                                    (w, j, v_b, i) in edges_total) - gp.quicksum(
            g_var[(v_b, i, w, j)] for w, j in vertices_total if
            (v_b, i, w, j) in edges_total) == 0)

    # Restriccion 5
    model.addConstrs(g_var[a, b, c, d] <= (len(neighborhoods) - 1) * y[a, b, c, d] for a, b, c, d in g_var.keys())
    # model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood] for v, i in vertices_barrier) == 1 for
    # vertices_neighborhood in vertices_neighborhood) model.addConstrs(gp.quicksum(y[vertices_neighborhood, v,
    # i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in vertices_neighborhood)
    #
    # model.addConstrs(gp.quicksum(g_var[v, i, vertices_neighborhood] for v, i in vertices_barrier) - gp.quicksum(g_var[
    # vertices_neighborhood, v, i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in
    # vertices_neighborhood) model.addConstrs(gp.quicksum(g_var[index] for index in edges_barrier if index[2] == v and
    # index[3] == i) + gp.quicksum(g_var[index] for index in edges_neighborhood if index[1] == v and index[2] == i) -
    # gp.quicksum(g_var[index] for index in edges_barrier if index[0] == v and index[1] == i) - gp.quicksum(g_var[index] for
    # index in edges_neighborhood if index[0] == v and index[1] == i) == 0 for v, i in vertices_barrier)

    # model.addConstrs(gp.quicksum(f[-1, w, k] for w in vertices_neighborhood if w <= -2) == 1 for k in
    # vertices_neighborhood if k <= -2) model.addConstrs(gp.quicksum(f[v, w, w] ))

    # model.addConstrs(gp.quicksum(z[w, v] for v in vertices_neighborhood if w != v) == 1 for w in
    # vertices_neighborhood)

    # flow conservation constraints for index in y_index: if len(index) == 3: model.addConstrs(gp.quicksum(y[tupla]
    # for tupla in edges_neighborhood if tupla[0] == v) == 1 for v in vertices_neighborhood)
    #
    # for v, i in vertices_barrier: tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[0] == v
    # and tupla[1] == i]) + gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v and tupla[1] ==
    # i]) tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[2] == v and tupla[3] == i]) +
    # gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[1] == v and tupla[2] == i])
    #
    #     model.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # for v in vertices_neighborhood:
    #     tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v])
    #     tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[2] == v])
    #
    #     model.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # model.addConstrs(gp.quicksum(y[tupla] for tupla in edges_neighborhood if tupla[2] == w) == 1 for w in
    # vertices_neighborhood)

    model.update()

    objective = gp.quicksum(dist[index] * y[index] for index in edges_total)

    model.setObjective(objective, GRB.MINIMIZE)

    second_time = time.time()

    time_elapsed = second_time - first_time

    model.update()

    model.Params.Threads = 6
    #        model.Params.timeLimit = time_limit - time_elapsed
    model.Params.timeLimit = time_limit
    model.Params.MIPGap = 0.05
    # model.Params.LazyConstraints = 1
    model.Params.NumericFocus = 1
    # model.Params.NonConvex = 2

    #        model.write('prueba.lp')
    #        model.write('prueba.mps')

    model.optimize()

    # model.computeIIS()
    # model.write('infactibilidad.ilp')
    # model.write('initial_sol.sol')

    results = [model.Runtime, model.ObjVal]

    y_indices = []

    for index in edges_total:
        if y[index].X > 0.5:
            y_indices.append(index)

    if log:
        print(y_indices)

    g_indices = {}

    for index in g_index:
        if g_var[index].X > 0.5:
            g_indices[index] = g_var[index].X

    if log:
        print(g_indices)

    # if picture:
    #     fig, ax = plt.subplots()
    #
    #     for b in barriers:
    #         ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')
    #
    #     for n in neighborhoods:
    #         ax.add_artist(n.artist)
    #
    #     p_vals = model.getAttr('x', point)
    #     print(p_vals)
    #
    #     points = []
    #     for keys, vals in p_vals.items():
    #         points.append(vals)
    #
    #     points = np.array(points).reshape((len(neighborhoods), 2))
    #     print(points)
    #
    #     for i in points:
    #         ax.scatter(i[0], i[1], s=10, c='black')
    #
    #     # print(points)
    #
    #     segments = []
    #
    #     for a, b, c, d in y_indices:
    #         if (a, b, c, d) in edges_neighborhood:
    #             if a < 0:
    #                 segments.append(
    #                     [points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
    #             if c < 0:
    #                 segments.append(
    #                     [barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
    #         if (a, b, c, d) in edges_barrier:
    #             segments.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])
    #
    #         # if (a, b, c, d) in ENN: segments.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][
    #         # 1], points[abs(c)-1][1]])
    #
    #     # print(segments)
    #     for segment in segments:
    #         ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
    #                  head_width=1, length_includes_head=True, color='black')
    #
    #     # plt.axis([-5, 105, -5, 105])
    #     plt.axis([0, 100, 0, 100])
    #
    #     ax.set_aspect('equal')
    #     plt.show()

    return results, y_indices, g_indices