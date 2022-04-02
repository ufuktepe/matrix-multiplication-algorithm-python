import sys
import math
import random
import time
import multiprocessing as mp

sys.setrecursionlimit(2000)

global THRESHOLD
global BASE_CASE

THRESHOLD = 30


def standard_multiply(a, b, c, n, m, idx):
    """
    Multiplies two matrices that are z-ordered.
    :param a: input matrix 1
    :param b: input matrix 2
    :param c: resultant matrix
    :param n: size of the input matrix
    :param m: size of the matrix for the base case
    :param idx: starting index for the resultant matrix
    :return: resultant matrix that is z-ordered
    """

    if n == m:
        for i in range(m):
            temp = [a[i * m + k] for k in range(m)]
            for j in range(m):
                r = 0
                for k in range(m):
                    r += temp[k] * b[j + k * m]
                c[idx] += r
                idx += 1

    else:
        p = len(a) // 4
        n = n // 2

        a11, a12, a21, a22 = partition(a)
        b11, b12, b21, b22 = partition(a)

        c = standard_multiply(a11, b11, c, n, m, idx)
        c = standard_multiply(a12, b21, c, n, m, idx)
        c = standard_multiply(a11, b12, c, n, m, idx + p)
        c = standard_multiply(a12, b22, c, n, m, idx + p)
        c = standard_multiply(a21, b11, c, n, m, idx + 2 * p)
        c = standard_multiply(a22, b21, c, n, m, idx + 2 * p)
        c = standard_multiply(a21, b12, c, n, m, idx + 3 * p)
        c = standard_multiply(a22, b22, c, n, m, idx + 3 * p)

    return c

def add(a, b):
    return [sum(i) for i in zip(a, b)]

def subtract(a, b):
    return [a_i - b_i for a_i, b_i in zip(a, b)]

def partition(a):
    """
    Partitions the given z-ordered matrix of size n x n into four n/2 x n/2 matrices.
    :param a: the matrix to be partitioned
    :return: partitioned matrices of size n/2 x n/2
    """
    p = len(a) // 4
    return a[:p], a[p:2*p], a[2*p:3*p], a[3*p:]


def strassen(a, b, n):
    if n <= THRESHOLD:
        c = [0] * len(a)
        return standard_multiply(a, b, c, n, m=BASE_CASE, idx=0)

    a11, a12, a21, a22 = partition(a)
    b11, b12, b21, b22 = partition(b)

    n = n // 2

    p1 = strassen(a11, subtract(b12, b22), n)
    p2 = strassen(add(a11, a12), b22, n)
    p3 = strassen(add(a21, a22), b11, n)
    p4 = strassen(a22, subtract(b21, b11), n)
    p5 = strassen(add(a11, a22), add(b11, b22), n)
    p6 = strassen(subtract(a12, a22), add(b21, b22), n)
    p7 = strassen(subtract(a11, a21), add(b11, b12), n)

    c11 = add(subtract(add(p5, p4), p2), p6)
    c12 = add(p1, p2)
    c21 = add(p3, p4)
    c22 = subtract(subtract(add(p5, p1), p3), p7)

    return c11 + c12 + c21 + c22


def pad(a):
    """
    Pads the given matrix with zeros. Repeatedly divides the dimension in half, each time taking the ceiling, until the
    dimension is less than or equal to the threshold. Then, repeatedly doubles the dimension until it hits or exceeds
    the original dimension in order to find the final dimension.
    :param a: matrix to be padded
    :return: the padded matrix and the size of the base case
    """

    n = len(a)

    if n <= THRESHOLD:
        return a, n

    # Repeatedly divide the dimension until it is less than or equal to the threshold
    m = n
    while m > THRESHOLD:
        m = math.ceil(m/2)

    base_case = m

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

    return a, base_case

def z_order(a_1D, a_2D, x, y, n, m):
    """
    Converts a 2D matrix (list of lists) into a z-ordered (morton-ordered) 1D list.
    :param a_1D: empty list to be populated (should be an empty list initially)
    :param a_2D: 2D matrix (list of lists) to be converted (input matrix)
    :param x: row index (should be 0 initially)
    :param y: column index (should be 0 initially)
    :param n: size of the 2D matrix
    :param m: size of the matrix for the base case
    :return: None
    """

    if n == m:
        for i in range(m):
            for j in range(m):
                a_1D.append(a_2D[y + i][x + j])
    else:
        n = n // 2
        z_order(a_1D, a_2D, x, y, n, m)
        z_order(a_1D, a_2D, x + n, y, n, m)
        z_order(a_1D, a_2D, x, y + n, n, m)
        z_order(a_1D, a_2D, x + n, y + n, n, m)


def reorder(c, d, n, m, row):
    """
    Converts a 1D z-ordered list into a 2D matrix.
    :param c: 1D z-ordered list
    :param d: resultant 2D matrix (list of lists)
    :param n: size of the matrix (square root of c)
    :param m: size of the matrix for the base case
    :param row: current row for the resultant 2D matrix (should be 0 initially)
    :return: resultant 2D matrix
    """

    if n == m:
        for i in range(m):
            for j in range(m):
                d[row + i].append(c[i * m + j])

    else:
        n = n // 2

        c11, c12, c21, c22 = partition(c)

        d = reorder(c11, d, n, m, row)
        d = reorder(c12, d, n, m, row)
        d = reorder(c21, d, n, m, row + n)
        d = reorder(c22, d, n, m, row + n)

    return d

def print_diagonal(a, n):
    for i in range(n):
        print(a[i][i])
    print('')

def read_input(d, input_file):
    a_2D = []
    b_2D = []
    row = []
    n_items = d ** 2

    count = 0
    with open(input_file) as f:
        for line in f:

            if line.strip() == '':
                continue
            count += 1

            try:
                row.append(int(float(line.strip())))
            except TypeError as e:
                print(e)

            if count % d == 0:
                if count <= n_items:
                    a_2D.append(row)
                else:
                    b_2D.append(row)
                row = []

    if len(a_2D) != d:
        raise ValueError('Invalid input!')
    elif len(b_2D) != d:
        raise ValueError('Invalid input!')

    return a_2D, b_2D


def run():

    d = int(sys.argv[2])
    input_file = sys.argv[3]

    a_2D, b_2D = read_input(d, input_file)

    global BASE_CASE

    a_2D, BASE_CASE = pad(a_2D)
    b_2D, BASE_CASE = pad(b_2D)

    # Size of the padded matrix
    padded_d = len(a_2D)

    # Create 1D arrays for Morton Ordering
    a_1D = []
    b_1D = []

    # Populate the 1D arrays
    z_order(a_1D, a_2D, 0, 0, padded_d, m=BASE_CASE)
    z_order(b_1D, b_2D, 0, 0, padded_d, m=BASE_CASE)

    # Multiply the matrices
    c = strassen(a_1D, b_1D, padded_d)

    # Reorder the resultant matrix as a 2D array
    c_ordered = [[] for i in range(padded_d)]
    c_ordered = reorder(c, c_ordered, padded_d, m=BASE_CASE, row=0)

    print_diagonal(c_ordered, d)


if __name__ == '__main__':
    run()

