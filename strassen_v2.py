import math


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
    a11 = [[a[i][j] for j in range(m)] for i in range(m)]

    # Top right matrix
    a12 = [[a[i][j] for j in range(m, n)] for i in range(m)]

    # Bottom left matrix
    a21 = [[a[i][j] for j in range(m)] for i in range(m, n)]

    # Bottom right matrix
    a22 = [[a[i][j] for j in range(m, n)] for i in range(m, n)]

    return a11, a12, a21, a22



def strassen(a, b):
    n = len(a)

    if n <= THRESHOLD:
        return standard_matrix_multiply(a, b)

    a11, a12, a21, a22 = partition(a)
    b11, b12, b21, b22 = partition(a)

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

    # Extend the existing rows with zeros
    for i in range(n):
        a[i].extend([0] * delta)

    # Add rows of zeros
    for i in range(delta):
        a.append([0] * m)

    return a




if __name__ == '__main__':

    THRESHOLD = 100

    zeros_large = [0] * 3

    a = [[1, 2, 3], [4, 5, 6]]

    for i in range(2):
        a.append(zeros_large)

    a[2][0] = 99

    prin