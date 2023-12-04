import numpy as np
import math

# DO NOT TOUCH THIS!
DIMENSION = 2

# READ ME!
# this program implements the newton interpolation algorithm for *any* polynomial
# interpolation to the power of n-1. --hugo

# TODO Make evaluation and find error for each p_n(x) where n = 1,2,3,...,n and x = 1,2,3,...

# PLEASE change the two constants below before running the program
N = 8  # this is the amount of (x,y) points
X_SUB = 2.5  # self-explanatory

# PLEASE define your function here
fx = lambda x: math.log(x)


def make_ys(n, fx, dim, start):
    # make function that wants to be interpolated
    x_test = np.zeros((n, dim))

    # smallest x value is 1
    for i in range(start, N + start):
        x_test[i - start, 0] = i
        # x_0 = 1
        x_test[i - start, 1] = fx(i)
        # x_test is a 2-col matrix
    return x_test


def newton_interpolation(n, fx, dim, start):
    # given n number of points of (x,y), find the n-1 amount of b coefficients iteratively

    # for fools who change constants that are not meant to be changed
    if dim != 2:
        print(
            "I changed back the value of DIMENSION to be = 2. Why? because you changed a constant that is not meant to be changed; I literally can not be bothered to make it be available for any dimensions nor should you care about it being applicable for any dimension."
        )
        dim = 2

    # init
    x_test = make_ys(n, fx, dim, start)
    b = np.zeros(n)

    # calculation columns init
    x_vals = np.arange(start, N + start, 1)
    col_now = np.zeros(n)
    col_prev = np.copy(x_test[:, 1])  # [0. 0.69314718 1.09861229 etc...]
    low = 0

    # assign f(x_0) to b_0
    b[0] = x_test[0, 1]

    for i in range(n - 1):
        # start from 1 to n-1 but i starts at 0
        col_now[:] = 0 #reset col_now
        k = 0
        for j in range(low, n - 1):
            col_now[j + 1] = (col_prev[j + 1] - col_prev[j]) / (
                x_vals[j + 1] - x_vals[k]
            )
            k += 1
        b[i + 1] = col_now[low + 1]
        col_prev = np.copy(col_now) #update col_prev with col_now
        # print(low)
        # print(col_now)
        low += 1

    return x_test, b


def evaluate(x, x_test, b, fx):
    # init
    n = len(b)
    x_true = x_test[:, 1]
    x_eval = np.zeros(n)
    high = 1

    # make array called x_eval where (x-1), (x-1)(x-2), (x-1)(x-2)(x-3), ...
    for i in range(n):
        product = 1
        for j in range(1, high):
            product *= x - j
            # print(j)
        x_eval[i] = product
        if high < n:
            high += 1

    # p_n = full-term interpolation eq.
    p_n = np.dot(b, x_eval)
    err = abs(fx(x) - p_n)

    return p_n, err

# TODO: cut-off function 

# get b coefficients
x_test, b = newton_interpolation(
    N, fx, DIMENSION, start=1
)  # why is the start value = 1? because log(0) = undefined.
print("b=", b)

# interpolate with (preferably) any non-integer x
p_n, err = evaluate(X_SUB, x_test, b, fx)
print(f"x = {X_SUB}\ny_actual = {fx(X_SUB)}\ny_interpolate = {p_n}\nerr = {err}")
