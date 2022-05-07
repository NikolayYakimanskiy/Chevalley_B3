#!/usr/bin/env python
# coding: utf-8

# In[85]:


# matrix of transformation from simple roots basis to standard R^3 basis
C = matrix([[1, -1, 0], [0, 1, -1], [0, 0, 1]])

# positive roots (first 3 are the simple roots)
# in simple roots basis
roots = [vector([1, 0, 0]), vector([0, 1, 0]), vector([0, 0, 1]), vector([1, 1, 0]), vector([0, 1, 1]),
         vector([2, 1, 0]), vector([1, 1, 1]), vector([2, 1, 1]), vector([2, 2, 1])]
# transforming to standard basis
for i in range(9):
    roots[i] = C * roots[i]
# negative roots
for i in range(9):
    roots.append(-roots[i])

dual_simple_roots = []
for i in range(3):
    a = roots[i]
    dual_simple_roots.append((2 / (a * a)) * a)
D = matrix(dual_simple_roots).transpose()


def is_root(root):
    return root in roots


def root_sign(root):
    if roots.index(root) < 9:
        return 1
    return -1


NN = [[0, -1, 0, -2, -1, 0, -2, 0, 0],
      [1, 0, -1, 0, 0, 0, 0, 1, 0],
      [0, 1, 0, 1, 0, 1, 0, 0, 0],
      [2, 0, -1, 0, 0, 0, -2, 0, 0],
      [1, 0, 0, 0, 0, -1, 0, 0, 0],
      [0, 0, -1, 0, 1, 0, 0, 0, 0],
      [2, 0, 0, 2, 0, 0, 0, 0, 0],
      [0, -1, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0]]

# N_ab = -N_ba = -N_-a,-b
# a + b + c = 0 ==> N_ab / (c,c) = N_bc / (a,a) = N_ca / (b,b)
# a + b + c + d = 0 and no pair is antipodal ==> N_ab*N_cd / (a+b,a+b) + N_bc*N_ad / (b+c,b+c) + N_ca*N_bd / (c+a,c+a) = 0
N = [[0] * 18 for _ in range(18)]
for i in range(18):
    for j in range(18):
        if not is_root(roots[i] + roots[j]):
            continue
        r = 0
        while is_root(roots[j] - r * roots[i]):
            r += 1
        N[i][j] = r

# determining signs...
# for s in range(9):
#    i0, j0 = -1, -1
#    for i in range(9):
#        if is_root(roots[s] - roots[i]) and root_sign(roots[s] - roots[i]) > 0:
#            # extraspecial pair for roots[s]
#            i0, j0 = i, roots.index(roots[s] - roots[i])
#            break
#    if i0 == -1:
#        continue
#
#    # special pairs
#    for i in range(i0 + 1, 9):
#        if is_root(roots[s] - roots[i]) and root_sign(roots[s] - roots[i]) > 0 and roots.index(roots[s] - roots[i]) > i:
#            j = roots.index(roots[s] - roots[i])
#            # roots[s] == roots[i0] + roots[j0] == roots[i] + roots[j]
#            t1 = 0
#            if is_root(roots[j] - roots[i0]):
#                k = roots.index(roots[j] - roots[i0])
#                t1 = ((roots[k] * roots[k]) / (roots[j] * roots[j])) * N[i0][k] * N[i][k]
#            t2 = 0
#            if is_root(roots[i] - roots[i0]):
#                k = roots.index(roots[i] - roots[i0])
#                t2 = ((roots[k] * roots[k]) / (roots[i] * roots[i])) * N[i0][k] * N[j][k]
#            N[i][j] *= sgn(t1 - t2)

# ...or just using precomputed ones
for i in range(9):
    for j in range(9):
        N[i][j] = NN[i][j]

# now signs are set for all pairs of roots a, b where 0 < a < b  (special pairs)
# setting signs for all pairs:
for i in range(18):
    for j in range(18):
        if not is_root(roots[i] + roots[j]):
            continue
        sg, ii, jj = 1, i, j
        if jj >= 9:
            ii, jj = roots.index(-roots[ii]), roots.index(-roots[jj])
            sg *= -1
        if ii >= 9:
            if root_sign(roots[ii] + roots[jj]) > 0:
                ii, jj = roots.index(-roots[ii]), roots.index(roots[ii] + roots[jj])
            else:
                ii, jj = jj, roots.index(-(roots[ii] + roots[jj]))
        # here roots[ii], roots[jj] > 0
        if ii > jj:
            ii, jj = jj, ii
            sg *= -1
        # roots[ii], roots[jj] is now a special pair (for which the sign of N is already set)
        N[i][j] = abs(N[i][j]) * sg * sgn(N[ii][jj])


def cartan(root1, root2):
    return 2 * (root1 * root2) / (root2 * root2)


def reflect(root_a, root_b):
    return root_b - cartan(root_b, root_a) * root_a


def ad(i, j):
    if i >= 18 and j >= 18:
        # ad[h_i, h_j] = 0
        return vector([0] * 21)
    if i >= 18 and j < 18:
        # ad[h_i, x_j] = <a_j, a_i> * x_j
        v = [0] * 21
        v[j] = cartan(roots[j], roots[i - 18])
        return vector(v)
    if i < 18 and j >= 18:
        # ad[x_i, h_j]
        return -ad(j, i)

    # i, j < 18
    if i == j:
        return vector([0] * 21)
    if roots[i] == -roots[j]:
        # ad[x_a, x_-a] = h_a
        a = roots[i]
        aa = (2 / (a * a)) * a  # dual root
        bb = D ^ (-1) * aa
        return vector([0] * 18 + list(bb))
    if is_root(roots[i] + roots[j]):
        v = [0] * 21
        v[roots.index(roots[i] + roots[j])] = N[i][j]
        return vector(v)
    else:
        return vector([0] * 21)


def ad_x(i):
    l = []
    for j in range(21):
        l.append(ad(i, j))
    return matrix(l).transpose()


X = []
for i in range(21):
    X.append(ad_x(i))


def ring_commute(a, b):
    return a * b - b * a


def group_commute(a, b):
    return a * b * a ^ (-1) * b ^ (-1)


def check_basis_relations():
    for i in range(18):
        for j in range(18):
            if is_root(roots[i] + roots[j]):
                assert ring_commute(X[i], X[j]) == N[i][j] * X[roots.index(roots[i] + roots[j])]
            elif roots[i] == -roots[j]:
                a = roots[i]
                aa = (2 / (a * a)) * a  # dual root
                bb = D ^ (-1) * aa
                assert ring_commute(X[i], X[j]) == bb[0] * X[18] + bb[1] * X[19] + bb[2] * X[20]
            else:
                assert ring_commute(X[i], X[j]) == 0
    for i in range(18, 21):
        for j in range(18, 21):
            assert ring_commute(X[i], X[j]) == 0
    for i in range(3):
        for j in range(18):
            assert ring_commute(X[18 + i], X[j]) == cartan(roots[j], roots[i]) * X[j]


def x(i, t=1):
    t_adx = t * ad_x(i)
    A = (t_adx) ^ 0
    M = A * 0
    k, fact_k = 0, 1
    while A != 0:
        M += A / fact_k
        A *= t_adx
        k += 1
        fact_k *= k
    return M


def w(i, t=1):
    # w_a(t) = x_a(t)x_-a(-t^-1)x_a(t)
    j = roots.index(-roots[i])
    x_i_t = x(i, t)
    return x_i_t * x(j, -(t ^ (-1))) * x_i_t


def h(i, t=1):
    # h_a(t) = w_a(t)w_a(1)^-1 = w_a(t)w_a(-1)
    return w(i, t) * w(i, -1)


def Q(i):
    # Q_a = w_a(1)x_a(1)
    return w(i) * x(i)


####################################################################################################

def matrix_of_variables(n, m, var_prefix):
    matr = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(var(var_prefix + "_" + str(i) + "_" + str(j)))
        matr.append(row)
    return matrix(matr)


def linearize_polynomial(p):
    q = p * 0  # not just 0 in order to get a polynomial object
    for mon, coeff in zip(p.monomials(), p.coefficients()):
        if mon.degree() < 2:
            q += coeff * mon
    return q


def inverse_and_linearize_const_var_matrix(a):
    # a = (a_const, a_var)
    # a_const is the constant part of matrix a (invertible!), a_var is the variable part; a = a_const + a_var
    a_const, a_var = a
    a_const_inv = a_const ^ (-1)
    # return a_const_inv * (1 - a_var * a_const_inv)
    return a_const_inv, -a_const_inv * a_var * a_const_inv


def multiply_and_linearize_const_var_matrices(a, b):
    # a = (a_const, a_var), b = (b_const, b_var)
    a_const, a_var = a
    b_const, b_var = b
    return a_const * b_const, a_const * b_var + a_var * b_const


def equations_mod2_from_matrix_equals_zero(a):
    # a is a symbolic matrix
    eqs = []
    for row in a:
        for elem in row:
            p = elem.polynomial(GF(2))
            assert p.constant_coefficient() == 0, "Polynomial has nonzero constant term"
            eqs.append(p)
    return eqs


def linearize_polynomials_and_construct_matrix(polys, vars2idx, start_row=0, mode='matrix'):
    # polys is a list of polynomials; vars2idx is a dict var_name (str) -> index (int)
    # allvars = wA.variables() + wB.variables() + wC.variables()
    # names2vars = {str(v): v for v in allvars}
    # vars2idx = {str(v): i for i, v in enumerate(allvars)}
    matr = {}  # to get a sparse matrix
    for i, p in enumerate(polys):
        assert p.constant_coefficient() == 0, "Polynomial has nonzero constant term"
        for mon in p.monomials():
            if mon.degree() == 1:
                matr[(start_row + i, vars2idx[str(mon)])] = GF(2)(1)
    if mode == 'matrix':
        return matrix(matr)
    else:  # mode == 'dict'
        return matr


####################################################################################################

def get_all_w_from_simples_and_linearize():
    st = set([0, 1, 2])
    ww = [None] * 18
    ww[0] = (w(0), matrix_of_variables(21, 21, 'A'))
    ww[1] = (w(1), matrix_of_variables(21, 21, 'B'))
    ww[2] = (w(2), matrix_of_variables(21, 21, 'C'))

    while len(st) < 18:
        for idx_a in st:
            for idx_b in st:
                new_idx = roots.index(reflect(roots[idx_a], roots[idx_b]))
                if ww[new_idx] is not None:
                    continue
                w_a, w_b, w_a_b = w(idx_a), w(idx_b), w(new_idx)
                ww[new_idx] = multiply_and_linearize_const_var_matrices(
                    multiply_and_linearize_const_var_matrices(ww[idx_a], ww[idx_b]),
                    inverse_and_linearize_const_var_matrix(ww[idx_a]))
                if w_a * w_b * w_a ^ (-1) == w_a_b:
                    pass
                elif w_a * w_b * w_a ^ (-1) == w_a_b ^ (-1):
                    ww[new_idx] = inverse_and_linearize_const_var_matrix(ww[new_idx])
                else:
                    assert False, "Something is wrong with relation R7 for w_{} and w_{}".format(idx_a, idx_b)
        for i in range(18):
            if ww[i] is not None:
                st.add(i)

    return ww


def get_all_x_from_simples_of_different_lengths_and_linearize(all_w, idxs=[0, 2]):
    # all_w is the result of get_all_w_from_simples()
    # idxs is the list of indices of two orthogonal simple roots of different lengths
    assert len(all_w) == 18
    assert len(idxs) == 2
    st = set(idxs)
    xx = [None] * 18
    for i in st:
        w_inv = inverse_and_linearize_const_var_matrix(all_w[i])
        xx[i] = (w_inv[0] * Q(i), w_inv[1] * Q(i))

    while len(st) < 18:
        for idx_a in range(18):
            for idx_b in st:
                new_idx = roots.index(reflect(roots[idx_a], roots[idx_b]))
                if xx[new_idx] is not None:
                    continue
                w_a, x_b, x_a_b = w(idx_a), x(idx_b), x(new_idx)
                xx[new_idx] = multiply_and_linearize_const_var_matrices(
                    multiply_and_linearize_const_var_matrices(all_w[idx_a], xx[idx_b]),
                    inverse_and_linearize_const_var_matrix(all_w[idx_a]))
                if w_a * x_b * w_a ^ (-1) == x_a_b:
                    pass
                elif w_a * x_b * w_a ^ (-1) == x_a_b ^ (-1):
                    xx[new_idx] = inverse_and_linearize_const_var_matrix(xx[new_idx])
                else:
                    assert False, "Something is wrong with relation R7 for w_{} and x_{}".format(idx_a, idx_b)
        for i in range(18):
            if xx[i] is not None:
                st.add(i)

    return xx


def get_R2_relations_params():
    # each r2_matr[a][b] is a list of triples (i, j, c_i_j) according to R2 relation for (x_a, x_b)
    r2_matr = [[None] * 18 for _ in range(18)]
    for a in range(18):
        for b in range(a, 18):
            if roots[a] + roots[b] == 0:
                continue
            lst = []
            for i in [1, 2]:
                for j in [1, 2]:
                    if is_root(i * roots[a] + j * roots[b]):
                        lst.append((i, j))
            assert len(lst) <= 2
            if len(lst) == 0:
                r2_matr[a][b] = []
            if len(lst) == 1:
                assert lst[0] == (1, 1)
                assert group_commute(x(a), x(b)) == x(roots.index(roots[a] + roots[b]), N[a][b]), "c_1_1 is not N_a_b"
                r2_matr[a][b] = [(1, 1, N[a][b])]
            if len(lst) == 2:
                assert lst[0] == (1, 1)
                i, j = lst[1]
                for c in [-1, 1]:
                    if group_commute(x(a), x(b)) == x(roots.index(roots[a] + roots[b]), N[a][b]) * x(
                            roots.index(i * roots[a] + j * roots[b]), c):
                        r2_matr[a][b] = [(1, 1, N[a][b]), (i, j, c)]
                        break
                assert r2_matr[a][b] is not None
    return r2_matr


def get_R2_R7_equations_linearize_mod2(all_w, all_x, r2_matr, vars2idx, limit=None, step=10):
    # all_w is the result of get_all_w_from_simples()
    # all_x is the result of get_all_x_from_simples_of_different_len()
    # r2_matr is the result of get_R2_relations_params()
    # vars2idx is a dict var_name (str) -> index (int)
    # limit is the "max" number of MATRIX equations (actually, some rough threshold after which we stop)
    # step is the step for recalculating the rank and adding corresponding equations if they increase it
    # the larger the step, the faster this function works and the more memory it requires!

    eqs = []
    matrix_eqs_cnt = 0
    max_rank = 0
    cur_eqs = []
    for a in range(18):
        for b in range(a + 1, 18):
            # ww (R7)
            new_idx = roots.index(reflect(roots[a], roots[b]))
            w_a, w_b, w_a_b = w(a), w(b), w(new_idx)
            cur_w = multiply_and_linearize_const_var_matrices(
                multiply_and_linearize_const_var_matrices(all_w[a], all_w[b]),
                inverse_and_linearize_const_var_matrix(all_w[a]))
            if w_a * w_b * w_a ^ (-1) == w_a_b:
                pass
            elif w_a * w_b * w_a ^ (-1) == w_a_b ^ (-1):
                cur_w = inverse_and_linearize_const_var_matrix(cur_w)
            else:
                assert False, "Something is wrong with relation R7 for w_{} and w_{}".format(a, b)
            matr = sum(cur_w) - sum(all_w[new_idx])
            if matr != 0:
                cur_eqs += equations_mod2_from_matrix_equals_zero(matr)

            # xx (R7)
            new_idx = roots.index(reflect(roots[a], roots[b]))
            w_a, x_b, x_a_b = w(a), x(b), x(new_idx)
            cur_x = multiply_and_linearize_const_var_matrices(
                multiply_and_linearize_const_var_matrices(all_w[a], all_x[b]),
                inverse_and_linearize_const_var_matrix(all_w[a]))
            if w_a * x_b * w_a ^ (-1) == x_a_b:
                pass
            elif w_a * x_b * w_a ^ (-1) == x_a_b ^ (-1):
                cur_x = inverse_and_linearize_const_var_matrix(cur_x)
            else:
                assert False, "Something is wrong with relation R7 for w_{} and x_{}".format(a, b)
            matr = sum(cur_x) - sum(all_x[new_idx])
            if matr != 0:
                cur_eqs += equations_mod2_from_matrix_equals_zero(matr)

            # R2
            if r2_matr[a][b] is not None:
                prod = matrix.identity(21)
                prod = (prod, 0 * prod)
                for i, j, c in r2_matr[a][b]:
                    assert c in [-2, -1, 1, 2]
                    tmp = all_x[roots.index(i * roots[a] + j * roots[b])]
                    if abs(c) > 1:  # abs(c) == 2
                        tmp = multiply_and_linearize_const_var_matrices(tmp, tmp)
                    if c < 0:
                        tmp = inverse_and_linearize_const_var_matrix(tmp)
                    prod = multiply_and_linearize_const_var_matrices(prod, tmp)
                xa_xb = multiply_and_linearize_const_var_matrices(all_x[a], all_x[b])
                xb_xa = multiply_and_linearize_const_var_matrices(all_x[b], all_x[a])
                matr = sum(xa_xb) - sum(multiply_and_linearize_const_var_matrices(prod, xb_xa))
                if matr != 0:
                    cur_eqs += equations_mod2_from_matrix_equals_zero(matr)

            matrix_eqs_cnt += 1
            if matrix_eqs_cnt % step == 0:
                rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
                if rk > max_rank:
                    eqs += cur_eqs
                    max_rank = rk
                cur_eqs = []

            if limit is not None and matrix_eqs_cnt > limit:
                rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
                if rk > max_rank:
                    eqs += cur_eqs
                    max_rank = rk
                return eqs

            print("matrix_eqs_cnt = {}, len = {}, max_rank = {}".format(matrix_eqs_cnt, len(eqs), max_rank))

    rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
    if rk > max_rank:
        eqs += cur_eqs
        max_rank = rk

    return eqs


####################################################################################################

def equations_mod2_from_matrix_in_var_t_equals_zero(a):
    # a is a symbolic matrix
    a = a.expand()
    eqs = []
    for row in a:
        for elem in row:
            p = elem.polynomial(GF(2))
            # ugly workaround because of some strange things in sage classes
            var_t = None
            for var in p.variables():
                if str(var) == 't':
                    var_t = var
                    break
            if var_t is None:
                assert p.constant_coefficient() == 0, "Polynomial has nonzero constant term"
                eqs.append(p)
            else:
                for pp in p.polynomial(var_t).coefficients():
                    # pp is a polynomial not containing t
                    assert pp.constant_coefficient() == 0, "Polynomial has nonzero constant term"
                    eqs.append(pp)
    return eqs


def get_all_xt_from_simples_of_different_lengths_and_linearize(t, idxs=[0, 2]):
    # t is ring element (variable)
    # idxs is the list of indices of two orthogonal simple roots of different lengths
    assert len(idxs) == 2
    st = set(idxs)
    x_t = [None] * 18
    for i, varname in zip(idxs, ['S', 'T']):
        x_t[i] = (x(i, t), matrix_of_variables(21, 21, varname))

    while len(st) < 18:
        for idx_a in range(18):
            for idx_b in st:
                new_idx = roots.index(reflect(roots[idx_a], roots[idx_b]))
                if x_t[new_idx] is not None:
                    continue
                w_a, x_b, x_a_b = w(idx_a), x(idx_b), x(new_idx)
                x_t[new_idx] = (w_a * x_t[idx_b][0] * w_a ^ (-1), w_a * x_t[idx_b][1] * w_a ^ (-1))
                if w_a * x_b * w_a ^ (-1) == x_a_b:
                    pass
                elif w_a * x_b * w_a ^ (-1) == x_a_b ^ (-1):
                    x_t[new_idx] = inverse_and_linearize_const_var_matrix(x_t[new_idx])
                else:
                    assert False, "Something is wrong with relation R7 for w_{} and x_{}".format(idx_a, idx_b)
        for i in range(18):
            if x_t[i] is not None:
                st.add(i)

    return x_t


def get_R2_R7_equations_for_xt_linearize_mod2(x_t, r2_matr, vars2idx, limit=None, step=10):
    # x_t is the result of get_all_xt_xtinv_from_simples_of_different_lengths_and_linearize()
    # r2_matr is the result of get_R2_relations_params()
    # vars2idx is a dict var_name (str) -> index (int)
    # limit is the "max" number of MATRIX equations (actually, some rough threshold after which we stop)
    # step is the step for recalculating the rank and adding corresponding equations if they increase it
    # the larger the step, the faster this function works and the more memory it requires!

    eqs = []
    matrix_eqs_cnt = 0
    max_rank = 0
    cur_eqs = []
    for a in range(18):
        for b in range(18):
            # x_t (R7)
            new_idx = roots.index(reflect(roots[a], roots[b]))
            w_a, x_b, x_a_b = w(a), x(b), x(new_idx)
            cur_x_t = (w_a * x_t[b][0] * w_a ^ (-1), w_a * x_t[b][1] * w_a ^ (-1))
            if w_a * x_b * w_a ^ (-1) == x_a_b:
                pass
            elif w_a * x_b * w_a ^ (-1) == x_a_b ^ (-1):
                cur_x_t = inverse_and_linearize_const_var_matrix(cur_x_t)
            else:
                assert False, "Something is wrong with relation R7 for w_{} and x_{}".format(a, b)
            matr = sum(cur_x_t) - sum(x_t[new_idx])
            if matr != 0:
                cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)

            # R2
            if r2_matr[a][b] is not None:
                x_a, x_b, x_a_t, x_b_t = (x(a), 0), (x(b), 0), x_t[a], x_t[b]

                if len(r2_matr[a][b]) == 0:
                    xa_xb = multiply_and_linearize_const_var_matrices(x_a, x_b_t)
                    xb_xa = multiply_and_linearize_const_var_matrices(x_b_t, x_a)
                    matr = sum(xa_xb) - sum(xb_xa)
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)

                    xa_xb = multiply_and_linearize_const_var_matrices(x_a_t, x_b)
                    xb_xa = multiply_and_linearize_const_var_matrices(x_b, x_a_t)
                    matr = sum(xa_xb) - sum(xb_xa)
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)

                    xa_xb = multiply_and_linearize_const_var_matrices(x_a_t, x_b_t)
                    xb_xa = multiply_and_linearize_const_var_matrices(x_b_t, x_a_t)
                    matr = sum(xa_xb) - sum(xb_xa)
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)
                elif len(r2_matr[a][b]) == 1:
                    i, j, c = r2_matr[a][b][0]
                    assert (i, j) == (1, 1)
                    prod = x_t[roots.index(i * roots[a] + j * roots[b])]
                    if abs(c) > 1:  # abs(c) == 2
                        prod = multiply_and_linearize_const_var_matrices(prod, prod)
                    if c < 0:
                        prod = inverse_and_linearize_const_var_matrix(prod)

                    xa_xb = multiply_and_linearize_const_var_matrices(x_a, x_b_t)
                    xb_xa = multiply_and_linearize_const_var_matrices(x_b_t, x_a)
                    matr = sum(xa_xb) - sum(multiply_and_linearize_const_var_matrices(prod, xb_xa))
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)

                    xa_xb = multiply_and_linearize_const_var_matrices(x_a_t, x_b)
                    xb_xa = multiply_and_linearize_const_var_matrices(x_b, x_a_t)
                    matr = sum(xa_xb) - sum(multiply_and_linearize_const_var_matrices(prod, xb_xa))
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)
                elif len(r2_matr[a][b]) == 2:
                    i, j, c = r2_matr[a][b][1]
                    assert (i, j) in [(1, 2), (2, 1)]
                    if i == 2:
                        xa_xb = multiply_and_linearize_const_var_matrices(x_a, x_b_t)
                        xb_xa = multiply_and_linearize_const_var_matrices(x_b_t, x_a)
                    elif j == 2:
                        xa_xb = multiply_and_linearize_const_var_matrices(x_a_t, x_b)
                        xb_xa = multiply_and_linearize_const_var_matrices(x_b, x_a_t)

                    prod = (matrix.identity(21), 0)
                    for i, j, c in r2_matr[a][b]:
                        assert c in [-2, -1, 1, 2]
                        tmp = x_t[roots.index(i * roots[a] + j * roots[b])]
                        if abs(c) > 1:  # abs(c) == 2
                            tmp = multiply_and_linearize_const_var_matrices(tmp, tmp)
                        if c < 0:
                            tmp = inverse_and_linearize_const_var_matrix(tmp)
                        prod = multiply_and_linearize_const_var_matrices(prod, tmp)
                    matr = sum(xa_xb) - sum(multiply_and_linearize_const_var_matrices(prod, xb_xa))
                    if matr != 0:
                        cur_eqs += equations_mod2_from_matrix_in_var_t_equals_zero(matr)

            matrix_eqs_cnt += 1
            if matrix_eqs_cnt % step == 0:
                rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
                if rk > max_rank:
                    eqs += cur_eqs
                    max_rank = rk
                cur_eqs = []

            if limit is not None and matrix_eqs_cnt > limit:
                rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
                if rk > max_rank:
                    eqs += cur_eqs
                    max_rank = rk
                return eqs

            print("matrix_eqs_cnt = {}, len = {}, max_rank = {}".format(matrix_eqs_cnt, len(eqs), max_rank))

    rk = linearize_polynomials_and_construct_matrix(eqs + cur_eqs, vars2idx).rank()
    if rk > max_rank:
        eqs += cur_eqs
        max_rank = rk

    return eqs


# In[78]:


ws = get_all_w_from_simples_and_linearize()
ws


# In[64]:


xs = get_all_x_from_simples_of_different_lengths_and_linearize(ws)
xs


# In[79]:


rs = get_R2_relations_params()
rs


# In[80]:


allvars = ws[0][1].variables() + ws[1][1].variables() + ws[2][1].variables()
vars2idx = {str(v): i for i, v in enumerate(allvars)}


# In[ ]:


vars2idx


# In[67]:


get_R2_R7_equations_linearize_mod2(ws, xs, rs, vars2idx)


# In[81]:


xt = get_all_xt_from_simples_of_different_lengths_and_linearize(1)
xt


# In[90]:


allvars2 = xt[0][1].variables() + xt[1][1].variables()
vars2idx2 = {str(v): i for i, v in enumerate(allvars2)}
vars2idx2


# In[ ]:


get_R2_R7_equations_for_xt_linearize_mod2(xt, rs, vars2idx2)


# In[ ]:




