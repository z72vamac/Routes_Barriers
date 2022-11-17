# TSPN-B

import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh

def sppn_b(barriers, neighborhoods, prepro=True, log=False, dominant=False, picture=False, time_limit=7200, init=False):
    if not (dominant):
        first_time = time.time()

        vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
        vertices_neighborhood = vertices_neighborhood[::-1]
        # list(itertools.product(vertices_neighborhood, range(len(barriers)), range(2))) + list(itertools.product(range(len(barriers)),
        # range(2), vertices_neighborhood))

        # for a, b in vertices_neighborhood:
        #     if type(neighborhoods[abs(a)-1]) is neigh.Circle:
        #         neighborhoods[abs(a)-1].radii = 1e-6
        #
        #     if type(neighborhoods[abs(a)-1]) is neigh.Poligonal:
        #         for dim in range(2):
        #             neighborhoods[abs(a)-1].V[1][dim] = neighborhoods[abs(a)-1].V[0][dim]

        edges_neighborhood = []

        for (a, b) in vertices_neighborhood:
            for c in range(len(barriers)):
                for d in range(2):
                    if prepro:
                        # Point of the barrier to check if is visible by the neighborhood
                        point = barriers[c][d]

                        # Neighborhood to check if it is visible by the point
                        neighborhood = neighborhoods[abs(a) - 1]

                        if af.cansee(point, neighborhood, barriers):
                            # Appending the feasible edges to edges_neighborhood
                            edges_neighborhood.append((a, b, c, d))
                            edges_neighborhood.append((c, d, a, b))
                    else:
                        edges_neighborhood.append((a, b, c, d))
                        edges_neighborhood.append((c, d, a, b))

        vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))

        edges_barrier = []
        for v, i in vertices_barrier:
            for w, j in vertices_barrier:
                if v < w:
                    # if prepro:
                    barrier = [barriers[v][i], barriers[w][j]]

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

        point_index = []
        for a, b in vertices_neighborhood:
            for dim in range(2):
                point_index.append((a, b, dim))

        alpha_index = []

        for a, b in vertices_barrier:
            for c, d, e, f in edges_neighborhood:
                if c < 0:
                    alpha_index.append((a, b, c, d, e, f))

        for a, b in vertices_total:
            for (c, d, e, f) in indices_barriers:
                alpha_index.append((a, b, c, d, e, f))

        if log:
            print("alpha = " + str(alpha_index))

        beta_index = []

        for a, b, c, d in indices_barriers:
            for e, f, g, h in edges_neighborhood:
                if e < 0:
                    beta_index.append((a, b, c, d, e, f, g, h))
                    beta_index.append((e, f, g, h, a, b, c, d))

        if log:
            print("beta = " + str(beta_index))

        gamma_index = beta_index

        delta_index = gamma_index

        epsilon_index = []

        for a, b, c, d in edges_neighborhood:
            if a < 0:
                epsilon_index.append((a, b, c, d))

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

        dist_index = edges_total

        if log:
            print("dist = " + str(dist_index))

        dif_index = []

        for a, b, c, d in edges_neighborhood:
            for dim in range(2):
                dif_index.append((a, b, c, d, dim))

        if log:
            print("dif = " + str(dif_index))

        # socp variables:
        d_inside_index = vertices_neighborhood

        dif_inside_index = []

        for (a, b) in vertices_neighborhood:
            for dim in range(2):
                dif_inside_index.append((a, b, dim))

        if log:
            print("d_inside_index = " + str(d_inside_index))
            print("dif_inside_index = " + str(dif_inside_index))

        # z variables:
        # z_index = list(itertools.product(vertices_neighborhood, vertices_neighborhood))
        # print("z_index = " + str(z_index))

        # f variables:
        g_index = y_index

        if log:
            print("g_index = " + str(g_index))

        model = gp.Model('HTSPS_Model')

        epsilon = model.addVars(epsilon_index, vtype=GRB.BINARY, name='epsilon')
        p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
        y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
        dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
        dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
        delta = model.addVars(delta_index, vtype=GRB.BINARY, name='delta')
        gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name='gamma')
        beta = model.addVars(beta_index, vtype=GRB.BINARY, name='beta')
        alpha = model.addVars(alpha_index, vtype=GRB.BINARY, name='alpha')

        # alpha00 = model.addVars(alpha00_index, vtype=GRB.BINARY, name='alpha00')
        # alpha10 = model.addVars(alpha10_index, vtype=GRB.BINARY, name='alpha10')
        # alpha_1 = model.addVars(alpha_1_index, vtype=GRB.BINARY, name='alpha_1')
        #
        # beta0 = model.addVars(beta0_index, vtype=GRB.BINARY, name='beta0')
        # beta1 = model.addVars(beta1_index, vtype=GRB.BINARY, name='beta1')

        # point = model.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'point')
        point = model.addVars(point_index, vtype=GRB.CONTINUOUS, name='point')

        d_inside = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='d_inside')
        dif_inside = model.addVars(dif_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif_inside')
        landa = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa')
        # z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
        g_var = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')

        model.update()

        # Alpha constraint
        L = -100000
        U = 100000

        # alpha-C
        for a, b, c, d, e, f in alpha_index:
            # Dos primeras
            # print((a, b, c, d, e, f))
            L = -100000
            U = 100000
            if (c, d, e, f) in indices_barriers:
                if a < 0:
                    if prepro:
                        barrera = [[barriers[c][d][0], barriers[c][d][1]], [barriers[e][f][0], barriers[e][f][1]]]
                        L, U = af.estima_det(neighborhoods[abs(a) - 1], barrera)

                        if U <= -1e-6:
                            alpha[a, b, c, d, e, f] = 0
                        elif L > 1e-6:
                            alpha[a, b, c, d, e, f] = 1
                        else:
                            model.addConstr(
                                (1 - alpha[a, b, c, d, e, f]) * L <= af.determinant([point[a, b, 0], point[a, b, 1]],
                                                                                    barriers[c][d], barriers[e][f]))
                            model.addConstr(
                                af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c][d], barriers[e][f]) <= U *
                                alpha[a, b, c, d, e, f])

                    else:
                        model.addConstr(
                            (1 - alpha[a, b, c, d, e, f]) * L <= af.determinant([point[a, b, 0], point[a, b, 1]],
                                                                                barriers[c][d], barriers[e][f]))
                        model.addConstr(
                            af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c][d], barriers[e][f]) <= U *
                            alpha[a, b, c, d, e, f])
                if a >= 0:
                    if prepro:
                        if af.determinant(barriers[a][b], barriers[c][d], barriers[e][f]) <= 0:
                            alpha[a, b, c, d, e, f] = 0
                        else:
                            alpha[a, b, c, d, e, f] = 1
                    else:
                        model.addConstr(
                            (1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b], barriers[c][d],
                                                                                barriers[e][f]))
                        model.addConstr(-U * alpha[a, b, c, d, e, f] <= -af.determinant(barriers[a][b], barriers[c][d],
                                                                                        barriers[e][f]))

            if a >= 0:
                if c < 0 and e >= 0:
                    if prepro:
                        barrera = [[barriers[a][b][0], barriers[a][b][1]], [barriers[e][f][0], barriers[e][f][1]]]
                        L, U = af.estima_det(neighborhoods[abs(c) - 1], barrera)

                        # print((L, U))

                        if U <= -1e-6:
                            alpha[a, b, c, d, e, f] = 0
                        elif L > 1e-6:
                            alpha[a, b, c, d, e, f] = 1
                        else:
                            L = -100000
                            U = 100000
                            # print('Hola')
                            # print((L, U))
                            model.addConstr((1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b],
                                                                                                [point[c, d, 0],
                                                                                                 point[c, d, 1]],
                                                                                                barriers[e][f]))
                            model.addConstr(
                                af.determinant(barriers[a][b], [point[c, d, 0], point[c, d, 1]], barriers[e][f]) <= U *
                                alpha[a, b, c, d, e, f])

                    else:
                        model.addConstr((1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b],
                                                                                            [point[c, d, 0],
                                                                                             point[c, d, 1]],
                                                                                            barriers[e][f]))
                        model.addConstr(
                            af.determinant(barriers[a][b], [point[c, d, 0], point[c, d, 1]], barriers[e][f]) <= U *
                            alpha[a, b, c, d, e, f])

                if c >= 0 and e < 0:
                    model.addConstr((1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b], barriers[c][d],
                                                                                        [point[e, f, 0],
                                                                                         point[e, f, 1]]))
                    model.addConstr(
                        af.determinant(barriers[a][b], barriers[c][d], [point[e, f, 0], point[e, f, 1]]) <= U * alpha[
                            a, b, c, d, e, f])
                if c < 0 and e < 0:
                    model.addConstr((1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b], [point[c, d, 0],
                                                                                                         point[
                                                                                                             c, d, 1]],
                                                                                        [point[e, f, 0],
                                                                                         point[e, f, 1]]))
                    model.addConstr(af.determinant(barriers[a][b], [point[c, d, 0], point[c, d, 1]],
                                                   [point[e, f, 0], point[e, f, 1]]) <= U * alpha[a, b, c, d, e, f])

        for a, b, c, d, e, f, g, h in beta_index:
            if ((a, b, c, d) in indices_barriers) or ((e, f, g, h) in indices_barriers):
                model.addConstr(
                    beta[a, b, c, d, e, f, g, h] == 2 * gamma[a, b, c, d, e, f, g, h] - alpha[a, b, e, f, g, h] - alpha[
                        c, d, e, f, g, h] + 1)

        for a, b, c, d, e, f, g, h in gamma_index:
            if ((a, b, c, d) in indices_barriers) or ((e, f, g, h) in indices_barriers):
                model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[a, b, e, f, g, h])
                model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[c, d, e, f, g, h])
                model.addConstr(gamma[a, b, c, d, e, f, g, h] >= alpha[a, b, e, f, g, h] + alpha[c, d, e, f, g, h] - 1)

        for a, b, c, d, e, f, g, h in delta_index:
            if (e, f, g, h) in indices_barriers:
                model.addConstr(0.5 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) <= delta[
                    a, b, c, d, e, f, g, h])
                model.addConstr(
                    2 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) >= delta[a, b, c, d, e, f, g, h])

        for a, b, c, d in epsilon_index:
            model.addConstr(
                gp.quicksum(delta[a, b, c, d, e, f, g, h] for e, f, g, h in indices_barriers) - len(barriers) + 1 <=
                epsilon[a, b, c, d])
            model.addConstr(len(barriers) * epsilon[a, b, c, d] <= gp.quicksum(
                delta[a, b, c, d, e, f, g, h] for e, f, g, h in indices_barriers))

            model.addConstr(y[a, b, c, d] + y[c, d, a, b] <= 2 * epsilon[a, b, c, d])

            # model.addConstr(epsilon[a, b, c, d] == epsilon[c, d, a, b])

            # if a < 0 and b < 0 and c == 0:
            #     model.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])

            # if no_ve(neighborhoods[abs(a)-1], barriers[b][0]) and no_ve(neighborhoods[abs(a)-1], barriers[b][1]):
            #     model.addConstr(epsilon[a, b, c] <= 0)

        # NS and NT constraints
        for a, b, dim in point.keys():
            neighborhood = neighborhoods[abs(a) - 1]

            if type(neighborhood) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - neighborhoods[abs(a) - 1].center[dim])
                model.addConstr(dif_inside[a, b, dim] >= neighborhoods[abs(a) - 1].center[dim] - point[a, b, dim])

                model.addConstr(
                    gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] *
                    d_inside[a, b])
                model.addConstr(d_inside[a, b] <= neighborhoods[abs(a) - 1].radii)
                # model.addConstr(point[a, b, dim] == neighborhoods[abs(a) - 1].center[dim])

            if type(neighborhood) is neigh.Poligonal:
                model.addConstrs(
                    point[a, b, dim] == landa[a, b] * neighborhood.V[0][dim] + (1 - landa[a, b]) *
                    neighborhood.V[1][dim] for dim in range(2))

        # model.addConstr(y[-1, 2, 1] >= 0.5)
        # model.addConstr(y[-2, 2, 1] >= 0.5)
        # model.addConstr(point[-1, 0] == 30)
        # model.addConstr(point[-1, 1] == 10)

        # model.addConstr(epsilon[-4, 3, 1] >= 0.5)
        # model.addConstr(epsilon[3, 1, -4] >= 0.5)

        # model.addConstr(y[3, 1, -1] >= 0.5)

        # dist constraints
        for a, b, c, d in dist_index:
            if (a, b, c, d) in edges_barrier:
                model.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d])))

            if (a, b, c, d) in edges_neighborhood:
                if a < 0:
                    model.addConstrs(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c][d][dim] for dim in range(2))
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c][d][dim] for dim in range(2))
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[
                            a, b, c, d] *
                        dist[a, b, c, d])
                if c < 0:
                    model.addConstrs(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a][b][dim] for dim in range(2))
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a][b][dim] for dim in range(2))
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[
                            a, b, c, d] *
                        dist[a, b, c, d])

            # if (a, b, c, d) in ENN: model.addConstr(dif[a, b, c, d, dim] >=   point[a, b, dim] - point[c, d,
            # dim]) model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + point[c, d, dim]) model.addConstr(
            # gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b,
            # c, d])

        l_out = 0
        u_out = 10000

        # p constraints
        for a, b, c, d in p.keys():

            if prepro:
                if a < 0:
                    neighborhood = neighborhoods[abs(a) - 1]
                    punto = barriers[c][d]
                    l_out = af.estima_L(neighborhood, punto)
                    u_out = af.estima_U(neighborhood, punto)

                if c < 0:
                    neighborhood = neighborhoods[abs(c) - 1]
                    punto = barriers[a][b]
                    l_out = af.estima_L(neighborhood, punto)
                    u_out = af.estima_U(neighborhood, punto)

            # print((l_out, u_out))
            model.addConstr(p[a, b, c, d] >= l_out * y[a, b, c, d])
            model.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - u_out * (1 - y[a, b, c, d]))

            # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)

        # model.addConstrs(z[v, v] == 0 for v in vertices_total

        for v, i in vertices_total:
            if v == -1:
                model.addConstr(gp.quicksum(y[v, i, w, j] for w, j in vertices_total if (v, i, w, j) in edges_total) - gp.quicksum(
                    y[w, j, v, i] for w, j in vertices_total if (w, j, v, i) in edges_total) == 1)
            elif v == -2:
                model.addConstr(gp.quicksum(y[v, i, w, j] for w, j in vertices_total if (v, i, w, j) in edges_total) - gp.quicksum(
                    y[w, j, v, i] for w, j in vertices_total if (w, j, v, i) in edges_total) == -1)

            else:
                model.addConstr(gp.quicksum(y[v, i, w, j] for w, j in vertices_total if (v, i, w, j) in edges_total) - gp.quicksum(
                    y[w, j, v, i] for w, j in vertices_total if (w, j, v, i) in edges_total) == 0)

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

        objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(
            dist[index] * y[index] for index in edges_barrier)

        model.setObjective(objective, GRB.MINIMIZE)

        second_time = time.time()

        time_elapsed = second_time - first_time

        model.update()

        model.Params.Threads = 6
        #        model.Params.timeLimit = time_limit - time_elapsed
        model.Params.timeLimit = time_limit

        # model.Params.LazyConstraints = 1
        model.Params.NumericFocus = 1
        # model.Params.NonConvex = 2

        #        model.write('prueba.lp')
        #        model.write('prueba.mps')

        model.optimize()

        results = [len(neighborhoods), len(barriers), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        if init:
            try:
                results[-2] = time_h
                results[-1] = objval_h
            except:
                print('No solution obtained by the heuristic')

        if model.Status == 3:
            model.computeIIS()
            model.write('infeasible_constraints.ilp')
            return results

        if model.SolCount == 0:
            return results

        #        model.write('solution.sol')

        results[2] = model.getAttr('MIPGap')
        results[3] = model.Runtime
        results[4] = time_elapsed
        results[5] = model.getAttr('NodeCount')
        results[6] = model.ObjVal

        y_indices = []

        for index in edges_total:
            if y[index].X > 0.5:
                y_indices.append(index)

        if log:
            print(y_indices)

        g_indices = []

        for index in g_index:
            if g_var[index].X > 0.5:
                g_indices.append(g_var[index])

        if log:
            print(g_indices)

        if picture:
            fig, ax = plt.subplots()

            for b in barriers:
                ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')

            for n in neighborhoods:
                ax.add_artist(n.artist)

            p_vals = model.getAttr('x', point)
            print(p_vals)

            points = []
            for keys, vals in p_vals.items():
                points.append(vals)

            points = np.array(points).reshape((len(neighborhoods), 2))
            print(points)

            for i in points:
                ax.scatter(i[0], i[1], s=10, c='black')

            # print(points)

            segments = []

            for a, b, c, d in y_indices:
                if (a, b, c, d) in edges_neighborhood:
                    if a < 0:
                        segments.append(
                            [points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
                    if c < 0:
                        segments.append(
                            [barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
                if (a, b, c, d) in edges_barrier:
                    segments.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])

                # if (a, b, c, d) in ENN: segments.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][
                # 1], points[abs(c)-1][1]])

            # print(segments)
            for segment in segments:
                ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
                         head_width=1, length_includes_head=True, color='black')

            # plt.axis([-5, 105, -5, 105])
            plt.axis([0, 100, 0, 100])

            ax.set_aspect('equal')
            plt.show()
