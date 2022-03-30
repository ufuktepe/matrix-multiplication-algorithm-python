import math
import random
import time


def standard_matrix_multiply(a, b):
    n = len(a)
    c = [[0 for y in range(n)] for x in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c[i][j] += a[i][k] * b[k][j]
    return c


def standard_matrix_multiply_caching(a, b):
    n = len(a)
    c = [[0 for y in range(n)] for x in range(n)]
    for i in range(n):
        for k in range(n):
            r = a[i][k]
            for j in range(n):
                c[i][j] += r * b[k][j]
    return c


def add(a, b):
    n = len(a)
    c = [[0 for y in range(n)] for x in range(n)]
    for i in range(n):
        for j in range(n):
            c[i][j] = a[i][j] + b[i][j]
    return c


def subtract(a, b):
    n = len(a)
    c = [[0 for y in range(n)] for x in range(n)]
    for i in range(n):
        for j in range(n):
            c[i][j] = a[i][j] - b[i][j]
    return c


def partition(a):
    """
    Partitions the given matrix of size n x n into four n/2 x n/2 matrices.
    :param a: the matrix to be partitioned
    :return: partitioned matrices of size n/2 x n/2
    """

    n = len(a)
    m = n // 2

    # Top left matrix
    # a11 = [[a_1D[i][j] for j in range(m)] for i in range(m)]
    # a11 = [a_1D[i][:m] for i in range(m)]

    a11 = []
    a12 = []
    a21 = []
    a22 = []

    for i in range(m):
        a11.append(a[i][:m])
        a12.append(a[i][m:n])

    for i in range(m, n):
        a21.append(a[i][:m])
        a22.append(a[i][m:n])


    # Top right matrix
    # a12 = [[a_1D[i][j] for j in range(m, n)] for i in range(m)]

    # a12 = [a_1D[i][m:n] for i in range(m)]

    # Bottom left matrix
    # a21 = [[a_1D[i][j] for j in range(m)] for i in range(m, n)]

    # a21 = [a_1D[i][:m] for i in range(m, n)]

    # Bottom right matrix
    # a22 = [[a_1D[i][j] for j in range(m, n)] for i in range(m, n)]

    # a22 = [a_1D[i][m:n] for i in range(m, n)]

    return a11, a12, a21, a22



def strassen(a, b):
    n = len(a)

    if n <= THRESHOLD:
        return standard_matrix_multiply_caching(a, b)

    a11, a12, a21, a22 = partition(a)
    b11, b12, b21, b22 = partition(b)

    p1 = strassen(a11, subtract(b12, b22))
    p2 = strassen(add(a11, a12), b22)
    p3 = strassen(add(a21, a22), b11)
    p4 = strassen(a22, subtract(b21, b11))
    p5 = strassen(add(a11, a22), add(b11, b22))
    p6 = strassen(subtract(a12, a22), add(b21, b22))
    p7 = strassen(subtract(a11, a21), add(b11, b12))

    c11 = add(subtract(add(p5, p4), p2), p6)
    c12 = add(p1, p2)
    c21 = add(p3, p4)
    c22 = subtract(subtract(add(p5, p1), p3), p7)

    # c12 = subtract(a21, a11)
    # c21 = add(b11, b12)
    # c22 = strassen(c12, c21)
    # c12 = subtract(a12, a22)
    # c21 = add(b21, b22)
    # c11 = strassen(c12, c21)
    # c12 = add(a11, a22)
    # c21 = add(b11, b22)
    # temp1 = strassen(c12, c21)
    # c11 = add(temp1, c11)
    # c22 = add(temp1, c22)
    # temp2 = add(a21, a22)
    # c21 = strassen(temp2, b11)
    # c22 = subtract(c22, c21)
    # temp1 = subtract(b21, b11)
    # temp2 = strassen(a22, temp1)
    # c21 = add(c21, temp2)
    # c11 = add(c11, temp2)
    # temp1 = subtract(b12, b22)
    # c12 = strassen(a11, temp1)
    # c22 = add(c22, c12)
    # temp2 = add(a11, a12)
    # temp1 = strassen(temp2, b22)
    # c12 = add(c12, temp1)
    # c11 = subtract(c11, temp1)

    # Construct the final matrix
    c = []
    for i in range(len(c11)):
        c.append(c11[i] + c12[i])
    for i in range(len(c21)):
        c.append(c21[i] + c22[i])

    return c


def pad(a):
    """
    Pads the given matrix with zeros. Repeatedly divides the dimension in half, each time taking the ceiling, until the
    dimension is less than or equal to the threshold. Then, repeatedly doubles the dimension until it hits or exceeds
    the original dimension in order to find the final dimension.
    :param a: matrix to be padded
    :return: the padded matrix
    """

    n = len(a)

    if n <= THRESHOLD:
        return a

    # Repeatedly divide the dimension until it is less than or equal to the threshold
    m = n
    while m > THRESHOLD:
        m = math.ceil(m/2)

    # Repeatedly double the dimension until it hits or exceeds the original dimension
    while m < n:
        m = m * 2

    # Calculate the number of extra rows/columns to be added
    delta = m - n

    if delta > 0:
        # Extend the existing rows with zeros
        for i in range(n):
            a[i].extend([0] * delta)

        # Add rows of zeros
        for i in range(delta):
            a.append([0] * m)

    return a


def run(d):
    a = []
    b = []

    for i in range(d):
        row = []
        for j in range(d):
            row.append(random.randint(0, 2))
        a.append(row)

    for i in range(d):
        row = []
        for j in range(d):
            row.append(random.randint(0, 2))
        b.append(row)

    # a = [[0, 1, 3, 2, 1, 2, 3, 1],
    #      [2, 1, 0, 3, 1, 2, 3, 1],
    #      [2, 3, 2, 1, 0, 2, 3, 1],
    #      [1, 0, 2, 2, 3, 1, 3, 1],
    #      [3, 1, 0, 1, 0, 1, 3, 1],
    #      [2, 2, 1, 1, 0, 3, 3, 1],
    #      [2, 1, 0, 3, 1, 2, 3, 1],
    #      [2, 3, 2, 1, 0, 2, 3, 1]]
    #
    # b = [[3, 1, 2, 1, 0, 3, 2, 0],
    #      [1, 0, 3, 3, 1, 2, 2, 0],
    #      [1, 3, 0, 1, 1, 3, 2, 0],
    #      [3, 2, 1, 2, 3, 2, 2, 0],
    #      [0, 3, 1, 2, 3, 0, 2, 0],
    #      [3, 1, 1, 3, 2, 2, 2, 0],
    #      [3, 1, 2, 1, 0, 3, 2, 0],
    #      [1, 0, 3, 3, 1, 2, 2, 0]]

    a = pad(a)
    b = pad(b)

    c = strassen(a, b)

    d = len(a)
    for i in range(d):
        for j in range(d):
            print(c[i][j], end=' ')
        print('')


if __name__ == '__main__':
    thresholds = []
    for i in range(16, 50):
        if i % 5 != 0:
            thresholds.append(i)

    dim = [2000]

    for d in dim:
        for t in thresholds:
            THRESHOLD = t
            start = time.time()
            run(d)
            end = time.time()
            print(f'{d},{t},{end - start}')

    # THRESHOLD = 100
    # start = time.time()
    # run(1600)
    # end = time.time()
    # print(f'{1600},{100},{end - start}')

