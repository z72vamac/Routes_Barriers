# TSPN-B

import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh
import estimacion_M as eM

from heuristic import heuristic
from heuristic2 import heuristic2


def tspn_b(
    barriers,
    neighborhoods,
    prepro=True,
    A4=True,
    dominant=False,
    log=False,
    picture=False,
    time_limit=7200,
    init=False,
):

    if not (dominant):
        first_time = time.time()

        # Indices of the vertices identified with the neighborhoods
        vertices_neighborhood = list(
            itertools.product(range(-len(neighborhoods), 0), range(1))
        )
        vertices_neighborhood = vertices_neighborhood[::-1]

        # Indices of the vertices identificed with the vertices of the barriers
        vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))

        # Indices of the edges joining a neighborhood and a barrier
        edges_neighborhood = []

        for (a, b) in vertices_neighborhood:
            for c, d in vertices_barrier:
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

        # Indices of the edges joining two barriers
        edges_barrier = []

        for v, i in vertices_barrier:
            for w, j in vertices_barrier:
                if v > w:
                    if prepro:
                        barrier = [barriers[v][i], barriers[w][j]]

                        intersect = False
                        for barrieri in barriers:
                            if af.intersect(barrieri, barrier):
                                intersect = True
                                break

                        if not (intersect):
                            edges_barrier.append((v, i, w, j))
                            edges_barrier.append((w, j, v, i))

                    else:
                        edges_barrier.append((v, i, w, j))
                        edges_barrier.append((w, j, v, i))

        indices_barriers = [(v, 0, v, 1) for v in range(len(barriers))]

        # Including edges joining two neighbourhoods
        edges_source_target = []

        if not (A4):
            for (a, b) in vertices_neighborhood:
                for (c, d) in vertices_neighborhood:
                    neighborhood1 = neighborhoods[abs(a) - 1]
                    neighborhood2 = neighborhoods[abs(c) - 1]

                    if af.canseeN(neighborhood1, neighborhood2, barriers):
                        edges_source_target.append((a, b, c, d))
                        # edges_source_target.append((c, d, a, b))

        # Completed graph
        vertices_total = vertices_neighborhood + vertices_barrier
        edges_total = edges_neighborhood + edges_barrier + edges_source_target

        if log:
            print("vertices_neighborhood = " + str(vertices_neighborhood))
            print("vertices_barrier = " + str(vertices_barrier))

            print("edges_neighborhood = " + str(edges_neighborhood))
            print("edges_barrier = " + str(edges_barrier))
            print("edges_source_target = " + str(edges_source_target))

        # Points in the neighborhoods
        point_index = []
        for a, b in vertices_neighborhood:
            for dim in range(2):
                point_index.append((a, b, dim))

        # Alpha index
        alpha_index = []

        for a, b in vertices_barrier:
            for c, d, e, f in edges_total:
                # if c < 0:
                alpha_index.append((a, b, c, d, e, f))

        for a, b in vertices_total:
            for (c, d, e, f) in indices_barriers:
                alpha_index.append((a, b, c, d, e, f))

        if log:
            print("alpha = " + str(alpha_index))

        beta_index = []

        for a, b, c, d in indices_barriers:
            for e, f, g, h in edges_total:
                beta_index.append((a, b, c, d, e, f, g, h))
                beta_index.append((e, f, g, h, a, b, c, d))

        if log:
            print("beta = " + str(beta_index))

        gamma_index = beta_index

        delta_index = gamma_index

        # delta_index = []
        #
        # for a, b, c, d in edges_neighborhood:
        #     if a < 0:
        #         for e, f, g, h in indices_barriers:
        #             delta_index.append((a, b, c, d, e, f, g, h))

        epsilon_index = []

        for a, b, c, d in edges_total:
            # if a < 0:
            epsilon_index.append((a, b, c, d))

        y_index = edges_total

        if log:
            print("y_index = " + str(y_index))

        p_index = edges_total

        if log:
            print("p_index = " + str(p_index))

        dist_index = edges_total

        if log:
            print("dist = " + str(dist_index))

        dif_index = []

        for a, b, c, d in edges_total:
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

        # Adjusting the bounds of the circles. Computing the biggest radii in the problem.
        # This is done to avoid the problem of the circles being too big.
        lb_min = 0

        lb_max = max(
            [ent.radii for ent in neighborhoods if type(ent) is neigh.Circle] + [lb_min]
        )

        model = gp.Model("Model: H-TSP-N")

        # Modeling the distance
        dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name="dist")
        dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name="dif")

        # Modeling determinants
        alpha = model.addVars(alpha_index, vtype=GRB.BINARY, name="alpha")
        beta = model.addVars(beta_index, vtype=GRB.BINARY, name="beta")
        gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name="gamma")
        delta = model.addVars(delta_index, vtype=GRB.BINARY, name="delta")
        epsilon = model.addVars(epsilon_index, vtype=GRB.BINARY, name="epsilon")

        # Modeling the conic neighborhoods
        point = model.addVars(point_index, vtype=GRB.CONTINUOUS, lb = -3*lb_max, name='point')
        d_inside = model.addVars(d_inside_index,
                                 vtype=GRB.CONTINUOUS,
                                 lb=0.0,
                                 name='d_inside')
        dif_inside = model.addVars(dif_inside_index,
                                   vtype=GRB.CONTINUOUS,
                                   lb=0.0,
                                   name='dif_inside')
        landa = model.addVars(d_inside_index,
                              vtype=GRB.CONTINUOUS,
                              lb=0.0,
                              ub=1.0,
                              name='landa')

        # Modeling the route
        y = model.addVars(y_index, vtype=GRB.BINARY, name="y")
        g_var = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name="g")

        # Modeling the product
        p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name="p")

        model.update()

        if init:
            # Creamos un heuristico como el de la k-mediana. Cogiendo los centros de los entornos.
            results_h, y_indices, g_indices = heuristic2(
                barriers,
                neighborhoods,
                prepro=True,
                log=False,
                dominant=False,
                picture=False,
            )

            time_h, objval_h = results_h[0], results_h[1]

            for a, b, c, d in y_indices:
                # model.addConstr(y[a, b, c, d] >= 0.5)
                y[a, b, c, d].start = 1

            for a, b, c, d in g_indices:
                g_var[a, b, c, d].start = g_indices[a, b, c, d]
            # model.read('initial_sol.sol')

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
                if (a, b) in vertices_neighborhood:
                    if prepro:
                        # Estimating L, U for this case
                        L, U = eM.estima_M_alpha1(
                            neighborhoods[abs(a) - 1], barriers[c][d], barriers[e][f]
                        )

                        # If U is negative, the whole neighborhood is in the hyperplane generated by the segment
                        if U < 0:
                            alpha[a, b, c, d, e, f] = 0

                        # If L is negative, the whole neighborhood is in the hyperplane generated by the segment
                        elif L > 0:
                            alpha[a, b, c, d, e, f] = 1
                        else:
                            model.addConstr(
                                (1 - alpha[a, b, c, d, e, f]) * L
                                <= af.determinant(
                                    [point[a, b, 0], point[a, b, 1]],
                                    barriers[c][d],
                                    barriers[e][f],
                                )
                            )
                            model.addConstr(
                                af.determinant(
                                    [point[a, b, 0], point[a, b, 1]],
                                    barriers[c][d],
                                    barriers[e][f],
                                )
                                <= U * alpha[a, b, c, d, e, f]
                            )
                    else:
                        model.addConstr(
                            (1 - alpha[a, b, c, d, e, f]) * L
                            <= af.determinant(
                                [point[a, b, 0], point[a, b, 1]],
                                barriers[c][d],
                                barriers[e][f],
                            )
                        )
                        model.addConstr(
                            af.determinant(
                                [point[a, b, 0], point[a, b, 1]],
                                barriers[c][d],
                                barriers[e][f],
                            )
                            <= U * alpha[a, b, c, d, e, f]
                        )

                elif (a, b) in vertices_barrier:
                    if (
                        af.determinant(barriers[a][b], barriers[c][d], barriers[e][f])
                        <= 0
                    ):
                        alpha[a, b, c, d, e, f] = 0
                    else:
                        alpha[a, b, c, d, e, f] = 1
                    # else:
                    #     model.addConstr((1 - alpha[a, b, c, d, e, f]) * L <= af.determinant(barriers[a][b], barriers[c][d], barriers[e][f]))
                    #     model.addConstr(-U * alpha[a, b, c, d, e, f] <= -af.determinant(barriers[a][b], barriers[c][d], barriers[e][f]))

            else:
                if (c, d, e, f) in edges_neighborhood:
                    if (c, d) in vertices_neighborhood:
                        if prepro:
                            L, U = eM.estima_M_alpha2(
                                barriers[a][b],
                                neighborhoods[abs(c) - 1],
                                barriers[e][f],
                            )

                            # If U is negative, the whole neighborhood is in the hyperplane generated by the segment
                            if U < 0:
                                alpha[a, b, c, d, e, f] = 0
                            # If L is negative, the whole neighborhood is in the hyperplane generated by the segment
                            elif L > 0:
                                alpha[a, b, c, d, e, f] = 1
                            else:
                                model.addConstr(
                                    (1 - alpha[a, b, c, d, e, f]) * L
                                    <= af.determinant(
                                        barriers[a][b],
                                        [point[c, d, 0], point[c, d, 1]],
                                        barriers[e][f],
                                    )
                                )
                                model.addConstr(
                                    af.determinant(
                                        barriers[a][b],
                                        [point[c, d, 0], point[c, d, 1]],
                                        barriers[e][f],
                                    )
                                    <= U * alpha[a, b, c, d, e, f]
                                )
                        else:
                            model.addConstr(
                                (1 - alpha[a, b, c, d, e, f]) * L
                                <= af.determinant(
                                    barriers[a][b],
                                    [point[c, d, 0], point[c, d, 1]],
                                    barriers[e][f],
                                )
                            )
                            model.addConstr(
                                af.determinant(
                                    barriers[a][b],
                                    [point[c, d, 0], point[c, d, 1]],
                                    barriers[e][f],
                                )
                                <= U * alpha[a, b, c, d, e, f]
                            )

                    elif (e, f) in vertices_neighborhood:
                        if prepro:
                            L, U = eM.estima_M_alpha3(
                                barriers[a][b],
                                barriers[c][d],
                                neighborhoods[abs(e) - 1],
                            )

                            # If U is negative, the whole neighborhood is in the hyperplane generated by the segment
                            if U < 0:
                                alpha[a, b, c, d, e, f] = 0

                            # If L is negative, the whole neighborhood is in the hyperplane generated by the segment
                            elif L > 0:
                                alpha[a, b, c, d, e, f] = 1
                            else:
                                model.addConstr(
                                    (1 - alpha[a, b, c, d, e, f]) * L
                                    <= af.determinant(
                                        barriers[a][b],
                                        barriers[c][d],
                                        [point[e, f, 0], point[e, f, 1]],
                                    )
                                )
                                model.addConstr(
                                    af.determinant(
                                        barriers[a][b],
                                        barriers[c][d],
                                        [point[e, f, 0], point[e, f, 1]],
                                    )
                                    <= U * alpha[a, b, c, d, e, f]
                                )
                        else:
                            model.addConstr(
                                (1 - alpha[a, b, c, d, e, f]) * L
                                <= af.determinant(
                                    barriers[a][b],
                                    barriers[c][d],
                                    [point[e, f, 0], point[e, f, 1]],
                                )
                            )
                            model.addConstr(
                                af.determinant(
                                    barriers[a][b],
                                    barriers[c][d],
                                    [point[e, f, 0], point[e, f, 1]],
                                )
                                <= U * alpha[a, b, c, d, e, f]
                            )

                elif (c, d, e, f) in edges_barrier:
                    if (
                        af.determinant(barriers[a][b], barriers[c][d], barriers[e][f])
                        <= 0
                    ):
                        alpha[a, b, c, d, e, f] = 0
                    else:
                        alpha[a, b, c, d, e, f] = 1

                elif (c, d, e, f) in edges_source_target:
                    if prepro:
                        L, U = eM.estima_M_alpha4(
                            barriers[a][b],
                            neighborhoods[abs(c) - 1],
                            neighborhoods[abs(e) - 1],
                        )

                        # If U is negative, the whole neighborhood is in the hyperplane generated by the segment
                        if U < 0:
                            alpha[a, b, c, d, e, f] = 0

                        # If L is negative, the whole neighborhood is in the hyperplane generated by the segment
                        elif L > 0:
                            alpha[a, b, c, d, e, f] = 1
                        else:
                            model.addConstr(
                                (1 - alpha[a, b, c, d, e, f]) * L
                                <= af.determinant(
                                    barriers[a][b],
                                    [point[c, d, 0], point[c, d, 1]],
                                    [point[e, f, 0], point[e, f, 1]],
                                )
                            )
                            model.addConstr(
                                af.determinant(
                                    barriers[a][b],
                                    [point[c, d, 0], point[c, d, 1]],
                                    [point[e, f, 0], point[e, f, 1]],
                                )
                                <= U * alpha[a, b, c, d, e, f]
                            )
                    else:
                        model.addConstr(
                            (1 - alpha[a, b, c, d, e, f]) * L
                            <= af.determinant(
                                barriers[a][b],
                                [point[c, d, 0], point[c, d, 1]],
                                [point[e, f, 0], point[e, f, 1]],
                            )
                        )
                        model.addConstr(
                            af.determinant(
                                barriers[a][b],
                                [point[c, d, 0], point[c, d, 1]],
                                [point[e, f, 0], point[e, f, 1]],
                            )
                            <= U * alpha[a, b, c, d, e, f]
                        )

        for a, b, c, d, e, f, g, h in beta_index:
            if (a, b, c, d) in indices_barriers or (e, f, g, h) in indices_barriers:
                model.addConstr(
                    beta[a, b, c, d, e, f, g, h]
                    == 2 * gamma[a, b, c, d, e, f, g, h]
                    - alpha[a, b, e, f, g, h]
                    - alpha[c, d, e, f, g, h]
                    + 1
                )
                model.addConstr(
                    gamma[a, b, c, d, e, f, g, h] <= alpha[a, b, e, f, g, h]
                )
                model.addConstr(
                    gamma[a, b, c, d, e, f, g, h] <= alpha[c, d, e, f, g, h]
                )
                model.addConstr(
                    gamma[a, b, c, d, e, f, g, h]
                    >= alpha[a, b, e, f, g, h] + alpha[c, d, e, f, g, h] - 1
                )

            if (e, f, g, h) in indices_barriers:

                model.addConstr(
                    0.5 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d])
                    <= delta[a, b, c, d, e, f, g, h]
                )
                model.addConstr(
                    2 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d])
                    >= delta[a, b, c, d, e, f, g, h]
                )

        for a, b, c, d in epsilon_index:
            model.addConstr(
                delta.sum(a, b, c, d, "*", "*", "*", "*") - len(barriers) + 1
                <= epsilon[a, b, c, d]
            )
            model.addConstr(
                len(barriers) * epsilon[a, b, c, d]
                <= delta.sum(a, b, c, d, "*", "*", "*", "*")
            )

            model.addConstr(y[a, b, c, d] + y[c, d, a, b] <= 2 * epsilon[a, b, c, d])

            # model.addConstr(len(barriers)*(y[a, b, c, d] + y[c, d, a, b]) <= 2 * gp.quicksum(delta[a, b, c, d, e, f, g, h] for e, f, g, h in indices_barriers))

            # model.addConstr(epsilon[a, b, c, d] == epsilon[c, d, a, b])

            # if a < 0 and b < 0 and c == 0:
            #     model.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])

            # if no_ve(neighborhoods[abs(a)-1], barriers[b][0]) and no_ve(neighborhoods[abs(a)-1], barriers[b][1]):
            #     model.addConstr(epsilon[a, b, c] <= 0)

        # N constraints
        for a, b, dim in point.keys():
            neighborhood = neighborhoods[abs(a) - 1]

            if type(neighborhood) is neigh.Circle:
                model.addConstr(
                    dif_inside[a, b, dim]
                    >= point[a, b, dim] - neighborhoods[abs(a) - 1].center[dim]
                )
                model.addConstr(
                    dif_inside[a, b, dim]
                    >= neighborhoods[abs(a) - 1].center[dim] - point[a, b, dim]
                )

                model.addConstr(
                    gp.quicksum(
                        dif_inside[a, b, dim] * dif_inside[a, b, dim]
                        for dim in range(2)
                    )
                    <= d_inside[a, b] * d_inside[a, b]
                )
                model.addConstr(d_inside[a, b] <= neighborhoods[abs(a) - 1].radii)

            if type(neighborhood) is neigh.Poligonal:
                model.addConstrs(
                    point[a, b, dim]
                    == landa[a, b] * neighborhood.V[0][dim]
                    + (1 - landa[a, b]) * neighborhood.V[1][dim]
                    for dim in range(2)
                )

        # dist constraints
        for a, b, c, d in dist_index:

            if (a, b, c, d) in edges_barrier:
                dist[a, b, c, d] = np.linalg.norm(
                    np.array(barriers[a][b]) - np.array(barriers[c][d])
                )

            elif (a, b, c, d) in edges_neighborhood:
                if (a, b) in vertices_neighborhood:
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c][d][dim]
                        for dim in range(2)
                    )
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= -point[a, b, dim] + barriers[c][d][dim]
                        for dim in range(2)
                    )
                    model.addConstr(
                        gp.quicksum(
                            dif[a, b, c, d, dim] * dif[a, b, c, d, dim]
                            for dim in range(2)
                        )
                        <= dist[a, b, c, d] * dist[a, b, c, d]
                    )

                elif (c, d) in vertices_neighborhood:
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a][b][dim]
                        for dim in range(2)
                    )
                    model.addConstrs(
                        dif[a, b, c, d, dim] >= -point[c, d, dim] + barriers[a][b][dim]
                        for dim in range(2)
                    )
                    model.addConstr(
                        gp.quicksum(
                            dif[a, b, c, d, dim] * dif[a, b, c, d, dim]
                            for dim in range(2)
                        )
                        <= dist[a, b, c, d] * dist[a, b, c, d]
                    )

            elif (a, b, c, d) in edges_source_target:
                model.addConstrs(
                    dif[a, b, c, d, dim] >= point[a, b, dim] - point[c, d, dim]
                    for dim in range(2)
                )
                model.addConstrs(
                    dif[a, b, c, d, dim] >= -point[a, b, dim] + point[c, d, dim]
                    for dim in range(2)
                )
                model.addConstr(
                    gp.quicksum(
                        dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)
                    )
                    <= dist[a, b, c, d] * dist[a, b, c, d]
                )

        l_out = 0
        u_out = 100000

        # p constraints
        for a, b, c, d in p.keys():

            if prepro:

                if (a, b, c, d) in edges_neighborhood:
                    if (a, b) in vertices_neighborhood:
                        neighborhood = neighborhoods[abs(a) - 1]
                        punto = neigh.Punto(barriers[c][d])
                        l_out, u_out = eM.estima_M_complete(neighborhood, punto)

                    if (c, d) in vertices_neighborhood:
                        neighborhood = neighborhoods[abs(c) - 1]
                        punto = neigh.Punto(barriers[a][b])
                        l_out, u_out = eM.estima_M_complete(neighborhood, punto)

                elif (a, b, c, d) in edges_barrier:
                    l_out = np.linalg.norm(
                        np.array(barriers[a][b]) - np.array(barriers[c][d])
                    )
                    u_out = np.linalg.norm(
                        np.array(barriers[a][b]) - np.array(barriers[c][d])
                    )

                elif (a, b, c, d) in edges_source_target:
                    neighborhood1 = neighborhoods[abs(a) - 1]
                    neighborhood2 = neighborhoods[abs(c) - 1]

                    l_out, u_out = eM.estima_M_complete(neighborhood1, neighborhood2)

            # print((l_out, u_out))
            model.addConstr(p[a, b, c, d] >= l_out * y[a, b, c, d])
            model.addConstr(
                p[a, b, c, d] >= dist[a, b, c, d] - u_out * (1 - y[a, b, c, d])
            )

            # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)

        # model.addConstrs(z[v, v] == 0 for v in vertices_total

        # Restriccion 1
        for v_n, j in vertices_neighborhood:
            model.addConstr(
                gp.quicksum(
                    y[v, i, v_n, j]
                    for v, i in vertices_total
                    if (v, i, v_n, j) in edges_neighborhood + edges_source_target
                )
                >= 1
            )

        # Restriccion 2
        for v, i in vertices_total:
            model.addConstr(
                gp.quicksum(
                    y[v, i, v_p, j]
                    for v_p, j in vertices_total
                    if (v, i, v_p, j) in edges_total
                )
                == gp.quicksum(
                    y[v_p, j, v, i]
                    for v_p, j in vertices_total
                    if (v_p, j, v, i) in edges_total
                )
            )

        # Restriccion 3
        for v_n, i in vertices_neighborhood:
            if v_n <= -2:
                model.addConstr(
                    gp.quicksum(
                        g_var[v_n, i, v, j]
                        for v, j in vertices_total
                        if (v_n, i, v, j) in edges_neighborhood + edges_source_target
                    )
                    - gp.quicksum(
                        g_var[v, j, v_n, i]
                        for v, j in vertices_total
                        if (v, j, v_n, i) in edges_neighborhood + edges_source_target
                    )
                    == 1
                )

        # Restriccion 4
        for v_b, i in vertices_barrier:
            model.addConstr(
                gp.quicksum(
                    g_var[(w, j, v_b, i)]
                    for w, j in vertices_total
                    if (w, j, v_b, i) in edges_total
                )
                - gp.quicksum(
                    g_var[(v_b, i, w, j)]
                    for w, j in vertices_total
                    if (v_b, i, w, j) in edges_total
                )
                == 0
            )

        # Restriccion 5
        model.addConstrs(
            g_var[a, b, c, d] <= (len(neighborhoods) - 1) * y[a, b, c, d]
            for a, b, c, d in g_var.keys()
        )
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

        objective = gp.quicksum(p[index] for index in p.keys())

        model.setObjective(objective, GRB.MINIMIZE)

        second_time = time.time()

        time_elapsed = second_time - first_time

        model.update()

        model.Params.Threads = 6
        #        model.Params.timeLimit = time_limit - time_elapsed
        model.Params.timeLimit = time_limit

        # model.Params.LazyConstraints = 1
        model.Params.NumericFocus = 1

        if not (A4):
            model.Params.NonConvex = 2

        # p = model.presolve()
        # p.write('presolved.lp')
        # model.write('prueba.lp')

        #        model.write('prueba.mps')

        model.optimize()

        results = [
            len(neighborhoods),
            len(barriers),
            prepro,
            A4,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]

        if init:
            try:
                results[9] = time_h
                results[10] = objval_h
            except:
                print("No solution obtained by the heuristic")

        if model.Status == 3:
            model.computeIIS()
            model.write("infeasible_constraints.ilp")
            return results

        if model.SolCount == 0:
            return results

        #        model.write('solution.sol')

        results[4] = model.getAttr("MIPGap")
        results[5] = model.Runtime
        results[6] = time_elapsed
        results[7] = model.getAttr("NodeCount")
        results[8] = model.ObjVal

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
                ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c="red")

            for n in neighborhoods:
                ax.add_artist(n.artist)

            p_vals = model.getAttr("x", point)
            print(p_vals)

            points = []
            for keys, vals in p_vals.items():
                points.append(vals)

            points = np.array(points).reshape((len(neighborhoods), 2))
            print(points)

            for i in points:
                ax.scatter(i[0], i[1], s=10, c="black")

            # print(points)

            segments = []

            for a, b, c, d in y_indices:
                if (a, b, c, d) in edges_neighborhood:
                    if (a, b) in vertices_neighborhood:
                        segments.append(
                            [
                                points[abs(a) - 1][0],
                                barriers[c][d][0],
                                points[abs(a) - 1][1],
                                barriers[c][d][1],
                            ]
                        )
                    if (c, d) in vertices_neighborhood:
                        segments.append(
                            [
                                barriers[a][b][0],
                                points[abs(c) - 1][0],
                                barriers[a][b][1],
                                points[abs(c) - 1][1],
                            ]
                        )
                if (a, b, c, d) in edges_barrier:
                    segments.append(
                        [
                            barriers[a][b][0],
                            barriers[c][d][0],
                            barriers[a][b][1],
                            barriers[c][d][1],
                        ]
                    )

                if (a, b, c, d) in edges_source_target:
                    segments.append(
                        [
                            points[abs(a) - 1][0],
                            points[abs(c) - 1][0],
                            points[abs(a) - 1][1],
                            points[abs(c) - 1][1],
                        ]
                    )
                # if (a, b, c, d) in ENN:

            # print(segments)
            for segment in segments:
                ax.arrow(
                    segment[0],
                    segment[2],
                    segment[1] - segment[0],
                    segment[3] - segment[2],
                    width=0.1,
                    head_width=1,
                    length_includes_head=True,
                    color="black",
                )

            # plt.axis([-5, 105, -5, 105])
            # plt.axis([0, 100, 0, 100])

            ax.set_aspect("equal")
            plt.show()
    # else:
    #     first_time = time.time()
    #
    #     vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
    #     vertices_neighborhood = vertices_neighborhood[::-1]
    #
    #     edges_neighborhood = []
    #
    #     for (a, b) in vertices_neighborhood:
    #         for c in range(len(barriers)):
    #             for d in range(2):
    #                 if prepro:
    #                     # Point of the barrier to check if is visible by the neighborhood
    #                     point = barriers[c][d]
    #
    #                     # Neighborhood to check if it is visible by the point
    #                     neighborhood = neighborhoods[abs(a) - 1]
    #
    #                     if af.cansee(point, neighborhood, barriers):
    #                         # Appending the feasible edges to edges_neighborhood
    #                         edges_neighborhood.append((a, b, c, d))
    #                         edges_neighborhood.append((c, d, a, b))
    #                 else:
    #                     edges_neighborhood.append((a, b, c, d))
    #                     edges_neighborhood.append((c, d, a, b))
    #
    #     vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))
    #
    #     edges_barrier = []
    #     for v, i in vertices_barrier:
    #         for w, j in vertices_barrier:
    #             if v != w:
    #                 if prepro:
    #                     barrier = [barriers[v][i], barriers[w][j]]
    #
    #                     if any([not (af.intersect(barrieri, barrier)) for barrieri in barriers]):
    #                         edges_barrier.append((v, i, w, j))
    #                 else:
    #                     edges_barrier.append((v, i, w, j))
    #
    #     vertices_total = vertices_neighborhood + vertices_barrier
    #     edges_total = edges_neighborhood + edges_barrier
    #
    #     if log:
    #         print("vertices_neighborhood = " + str(vertices_neighborhood))
    #         print("vertices_barrier = " + str(vertices_barrier))
    #
    #         print("edges_neighborhood = " + str(edges_neighborhood))
    #         print("edges_barrier = " + str(edges_barrier))
    #
    #     point_index = []
    #     for a, b in vertices_neighborhood:
    #         for dim in range(2):
    #             point_index.append((a, b, dim))
    #
    #     y_index = edges_total
    #
    #     if log:
    #         print("y = " + str(y_index))
    #
    #     dist_index = edges_total
    #
    #     if log:
    #         print("dist = " + str(dist_index))
    #
    #     dif_index = []
    #
    #     for a, b, c, d in edges_total:
    #         for dim in range(2):
    #             dif_index.append((a, b, c, d, dim))
    #
    #     # P_S and P_T: indices of the points in the neighborhoods
    #     p_index = edges_neighborhood
    #
    #     if log:
    #         print("p_index = " + str(p_index))
    #
    #     dom_set = af.dominant_set(neighborhoods, barriers)
    #
    #     # z variables:
    #     z_index = list(dom_set.keys())
    #     print("z_index = " + str(z_index))
    #
    #     # f variables:
    #     g_index = y_index
    #
    #     if log:
    #         print("g_index = " + str(g_index))
    #
    #
    #     model = gp.Model('HTSPN_Model_dominant')
    #
    #     p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
    #     y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
    #     dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
    #     dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
    #
    #     # point = model.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'point')
    #     point = model.addVars(point_index, vtype=GRB.CONTINUOUS, name='point')
    #
    #     z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    #     g_var = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g_var')
    #
    #     model.update()
    #
    #     if init:
    #         time_h, objval_h = heuristic(barriers, neighborhoods)
    #
    #     # NS and NT constraints
    #     for a, b, dim in point.keys():
    #
    #         model.addConstr(point[a, b, dim] == gp.quicksum(dom_set[index][dim]*z[index] for index in dom_set.keys() if index[0] == a))
    #
    #
    #     for a, c, d, e, f in z.keys():
    #         model.addConstr(z[a, c, d, e, f] <= y[c, d, a, 0])
    #         model.addConstr(z[a, c, d, e, f] <= y[a, 0, e, f])
    #         model.addConstr(z[a, c, d, e, f] >= y[c, d, a, 0] + y[a, 0, e, f] - 1)
    #
    #     # z constraints
    #     # for
    #     # model.addConstr(y[-1, 2, 1] >= 0.5)
    #     # model.addConstr(y[-2, 2, 1] >= 0.5)
    #     # model.addConstr(point[-1, 0] == 30)
    #     # model.addConstr(point[-1, 1] == 10)
    #
    #     # model.addConstr(epsilon[-4, 3, 1] >= 0.5)
    #     # model.addConstr(epsilon[3, 1, -4] >= 0.5)
    #
    #     # model.addConstr(y[3, 1, -1] >= 0.5)
    #
    #     # dist constraints
    #     for a, b, c, d, dim in dif_index:
    #         if (a, b, c, d) in edges_barrier:
    #             model.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d])))
    #
    #         if (a, b, c, d) in edges_neighborhood:
    #             if a < 0:
    #                 model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c][d][dim])
    #                 model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c][d][dim])
    #                 model.addConstr(
    #                     gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
    #                     dist[a, b, c, d])
    #             if c < 0:
    #                 model.addConstr(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a][b][dim])
    #                 model.addConstr(dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a][b][dim])
    #                 model.addConstr(
    #                     gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
    #                     dist[a, b, c, d])
    #
    #         # if (a, b, c, d) in ENN: model.addConstr(dif[a, b, c, d, dim] >=   point[a, b, dim] - point[c, d,
    #         # dim]) model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + point[c, d, dim]) model.addConstr(
    #         # gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b,
    #         # c, d])
    #
    #     l_out = 0
    #     u_out = 10000
    #
    #     # p constraints
    #     for a, b, c, d in p.keys():
    #
    #         if a < 0:
    #             neighborhood = neighborhoods[abs(a) - 1]
    #             punto = barriers[c][d]
    #             l_out = af.estima_L(neighborhood, punto)
    #             u_out = af.estima_U(neighborhood, punto)
    #
    #         if c < 0:
    #             neighborhood = neighborhoods[abs(c) - 1]
    #             punto = barriers[a][b]
    #             l_out = af.estima_L(neighborhood, punto)
    #             u_out = af.estima_U(neighborhood, punto)
    #
    #         # print((l_out, u_out))
    #         model.addConstr(p[a, b, c, d] >= l_out * y[a, b, c, d])
    #         model.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - u_out * (1 - y[a, b, c, d]))
    #
    #         # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)
    #
    #     # model.addConstrs(z[v, v] == 0 for v in vertices_total
    #
    #     # Restriccion 1
    #     model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood, j] for v, i in vertices_total if
    #                                  (v, i, vertices_neighborhood, j) in edges_neighborhood) >= 1 for
    #                      vertices_neighborhood, j in vertices_neighborhood)
    #
    #     # Restriccion 2
    #     for v, i in vertices_total:
    #         model.addConstr(gp.quicksum(y[v, i, vertices_neighborhood, j] for vertices_neighborhood, j in vertices_total if
    #                                     (v, i, vertices_neighborhood, j) in edges_total) == gp.quicksum(
    #             y[vertices_neighborhood, j, v, i] for vertices_neighborhood, j in vertices_total if
    #             (vertices_neighborhood, j, v, i) in edges_total))
    #
    #     # Restriccion 3
    #     for vertices_neighborhood, i in vertices_neighborhood:
    #         if vertices_neighborhood <= -2:
    #             model.addConstr(gp.quicksum(g_var[vertices_neighborhood, i, v, j] for v, j in vertices_total if
    #                                         (vertices_neighborhood, i, v, j) in edges_neighborhood) - gp.quicksum(
    #                 g_var[v, j, vertices_neighborhood, i] for v, j in vertices_total if
    #                 (v, j, vertices_neighborhood, i) in edges_neighborhood) == 1)
    #
    #     # Restriccion 4
    #     for vertices_barrier, i in vertices_barrier:
    #         model.addConstr(gp.quicksum(g_var[(w, j, vertices_barrier, i)] for w, j in vertices_total if
    #                                     (w, j, vertices_barrier, i) in edges_total) - gp.quicksum(
    #             g_var[(vertices_barrier, i, w, j)] for w, j in vertices_total if
    #             (vertices_barrier, i, w, j) in edges_total) == 0)
    #
    #     # Restriccion 5
    #     model.addConstrs(g_var[a, b, c, d] <= (len(neighborhoods) - 1) * y[a, b, c, d] for a, b, c, d in g_var.keys())
    #     # model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood] for v, i in vertices_barrier) == 1 for
    #     # vertices_neighborhood in vertices_neighborhood) model.addConstrs(gp.quicksum(y[vertices_neighborhood, v,
    #     # i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in vertices_neighborhood)
    #     #
    #     # model.addConstrs(gp.quicksum(g_var[v, i, vertices_neighborhood] for v, i in vertices_barrier) - gp.quicksum(g_var[
    #     # vertices_neighborhood, v, i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in
    #     # vertices_neighborhood) model.addConstrs(gp.quicksum(g_var[index] for index in edges_barrier if index[2] == v and
    #     # index[3] == i) + gp.quicksum(g_var[index] for index in edges_neighborhood if index[1] == v and index[2] == i) -
    #     # gp.quicksum(g_var[index] for index in edges_barrier if index[0] == v and index[1] == i) - gp.quicksum(g_var[index] for
    #     # index in edges_neighborhood if index[0] == v and index[1] == i) == 0 for v, i in vertices_barrier)
    #
    #     # model.addConstrs(gp.quicksum(f[-1, w, k] for w in vertices_neighborhood if w <= -2) == 1 for k in
    #     # vertices_neighborhood if k <= -2) model.addConstrs(gp.quicksum(f[v, w, w] ))
    #
    #     # model.addConstrs(gp.quicksum(z[w, v] for v in vertices_neighborhood if w != v) == 1 for w in
    #     # vertices_neighborhood)
    #
    #     # flow conservation constraints for index in y_index: if len(index) == 3: model.addConstrs(gp.quicksum(y[tupla]
    #     # for tupla in edges_neighborhood if tupla[0] == v) == 1 for v in vertices_neighborhood)
    #     #
    #     # for v, i in vertices_barrier: tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[0] == v
    #     # and tupla[1] == i]) + gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v and tupla[1] ==
    #     # i]) tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[2] == v and tupla[3] == i]) +
    #     # gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[1] == v and tupla[2] == i])
    #     #
    #     #     model.addConstr(tuplas_salen - tuplas_entran == 0)
    #     #
    #     # for v in vertices_neighborhood:
    #     #     tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v])
    #     #     tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[2] == v])
    #     #
    #     #     model.addConstr(tuplas_salen - tuplas_entran == 0)
    #     #
    #     # model.addConstrs(gp.quicksum(y[tupla] for tupla in edges_neighborhood if tupla[2] == w) == 1 for w in
    #     # vertices_neighborhood)
    #
    #     model.update()
    #
    #     objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(
    #         dist[index] * y[index] for index in edges_barrier)
    #     model.setObjective(objective, GRB.MINIMIZE)
    #
    #     second_time = time.time()
    #
    #     time_elapsed = second_time - first_time
    #
    #     model.update()
    #
    #     model.Params.Threads = 6
    #     model.Params.timeLimit = time_limit - time_elapsed
    #     # model.Params.LazyConstraints = 1
    #     model.Params.NumericFocus = 1
    #     # model.Params.NonConvex = 2
    #
    #     model.write('prueba.lp')
    #     model.write('prueba.mps')
    #
    #     model.optimize()
    #
    #     results = [len(neighborhoods), len(barriers), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    #
    #     if init:
    #         try:
    #             results[-2] = time_h
    #             results[-1] = objval_h
    #         except:
    #             print('No solution obtained by the heuristic')
    #
    #
    #
    #     if model.SolCount == 0:
    #         return results
    #
    #     model.write('solution.sol')
    #
    #     results[2] = model.getAttr('MIPGap')
    #     results[3] = model.Runtime + time_elapsed
    #     results[4] = model.getAttr('NodeCount')
    #     results[5] = model.ObjVal
    #
    #     y_indices = []
    #
    #     for index in edges_total:
    #         if y[index].X > 0.5:
    #             y_indices.append(index)
    #
    #     if log:
    #         print(y_indices)
    #
    #     g_indices = []
    #
    #     for index in g_index:
    #         if g_var[index].X > 0.5:
    #             g_indices.append(g_var[index])
    #
    #     if log:
    #         print(g_indices)
    #
    #     if picture:
    #         fig, ax = plt.subplots()
    #
    #         for b in barriers:
    #             ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')
    #
    #         for n in neighborhoods:
    #             ax.add_artist(n.artist)
    #
    #         p_vals = model.getAttr('x', point)
    #         print(p_vals)
    #
    #         points = []
    #         for keys, vals in p_vals.items():
    #             points.append(vals)
    #
    #         points = np.array(points).reshape((len(neighborhoods), 2))
    #         print(points)
    #
    #         for i in points:
    #             ax.scatter(i[0], i[1], s=10, c='black')
    #
    #         # print(points)
    #
    #         segments = []
    #
    #         for a, b, c, d in y_indices:
    #             if (a, b, c, d) in edges_neighborhood:
    #                 if a < 0:
    #                     segments.append(
    #                         [points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
    #                 if c < 0:
    #                     segments.append(
    #                         [barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
    #             if (a, b, c, d) in edges_barrier:
    #                 segments.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])
    #
    #             # if (a, b, c, d) in ENN: segments.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][
    #             # 1], points[abs(c)-1][1]])
    #
    #         # print(segments)
    #         for segment in segments:
    #             ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
    #                      head_width=1, length_includes_head=True, color='black')
    #
    #         # plt.axis([-5, 105, -5, 105])
    #         plt.axis([0, 100, 0, 100])
    #
    #         ax.set_aspect('equal')
    #         plt.show()
    #
    #     pass

    print(results)
    return results
