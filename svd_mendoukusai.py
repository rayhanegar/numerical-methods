import numpy as np

# example question:
# a = np.array([[3, 1, 1], [-1, 3, 1]])
# a = np.array([[3, 2, 2], [2, 3, -2]])
# a = np.array([[3, 2, 0], [2, 0, 0]])

a = np.array([[2, 4], [1, 3], [0, 0], [0, 0]])

# matrix A of size mxn where A = U * S * V.T
# U is an orthogonal matrix of size mxn
# S is a diagonal matrix of size mxn is the 2nd power of
# singular eigenvalues whereby if you square root S you will
# get the singular eigenvalues
# V^T is an orthogonal matrix of size nxn is the eigenvectors

# A.T * A = V * S * S.T * V.T
# A * A.T = U * S * S.T * U.T
# A * V = V * S^2

# implementation of svd is not iterative

def norm(vec):
    norm = 0
    for i in range(len(vec)):
        norm += vec[i] ** 2
    return norm

def cur_vec(a, i):
    vec = a[0 : len(a), i]
    return vec

def gs2(a):
    # initialise final e vector
    SIZE = len(a)
    e = np.zeros((SIZE, SIZE))

    # k = 1
    # grab current vector
    vec = cur_vec(a, 0)

    # calculate norm
    temp = norm(vec)

    # find e_1
    for i in range(SIZE):
        e[i, 0] = vec[i] / np.sqrt(temp)

    # k > 1
    for k in range(1, SIZE):
        vec = cur_vec(a, k)
        temp_res = 0

        for l in range(0, k):
            vec_e = cur_vec(e, l)
            temp_res += np.dot(vec, vec_e) * vec_e
        numerator = vec - temp_res

        for m in range(SIZE):
            e[m, k] = numerator[m] / np.sqrt(norm(numerator))

    return e

def qr(a, iteration=16):
    # this function uses QR decomposition
    a = np.copy(a)
    n = len(a)
    val = np.zeros((n, n))

    # inject 1 to columns without a leading 1
    for k in range(n):
        cvc = a[:, k]
        sum = 0
        for i in range(len(cvc)):
            sum += cvc[i]
        if sum > 0:
            val[k] = cvc
        elif sum == 0:
            cvc[k] = 1
            val[k] = cvc
    val = val.T
    # print("val=", val)

    vec = np.eye(n)

    for _ in range(iteration):
        Q = gs2(val)
        vec = vec @ np.copy(Q)
        R = Q.T @ val  # A_k = Q_k * R_k
        val = R @ Q  # A_(k+1) = R_k * Q_k

    return val.diagonal(), vec


def svd(a):
    A = np.copy(a)

    SL = A @ A.T
    SR = A.T @ A

    val_L, vec_L = qr(SL)
    val_R, vec_R = qr(SR)

    U = np.zeros((len(SL), len(SL)))
    for k in range(len(SL)):
        cvc = vec_L[:, k]

        norm = 0
        for i in range(len(cvc)):
            norm += cvc[i] ** 2
        norm = np.sqrt(norm)

        ncvc = np.zeros(len(SL))
        for i in range(len(SL)):
            ncvc[i] = cvc[i] / norm

        U[k] = ncvc
    U = U.T 

    V = np.zeros((len(SR), len(SR)))
    for k in range(len(SR)):
        cvc = vec_R[:, k]

        norm = 0
        for i in range(len(cvc)):
            norm += cvc[i] ** 2
        norm = np.sqrt(norm)

        ncvc = np.zeros(len(SR))
        for i in range(len(SR)):
            ncvc[i] = cvc[i] / norm

        V[k] = ncvc

    size = np.shape(A)
    S = np.zeros(size)

    for i in range(min(size[0], size[1])):
        S[i, i] = np.sqrt(val_L[i])

    return U, S, V

print(">>> soal")
print("a=", a)

print("\n>>> numpy.linalg.svd")
u, s, v = np.linalg.svd(a)
print("u=\n", np.round(u, decimals=5))
print("s=\n", np.round(s, decimals=5))
print("v=\n", np.round(v, decimals=5))

print("\n>>> milik kelompok 3")
u2, s2, v2 = svd(a)
print("u=\n", np.round(u2, decimals=5))
print("s=\n", np.round(s2, decimals=5))
print("v=\n", np.round(v2, decimals=5))
