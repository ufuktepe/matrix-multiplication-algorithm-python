import sys
import math
import random
import time
import multiprocessing as mp

sys.setrecursionlimit(2000)

global THRESHOLD
global BASE_CASE


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

def print_matrix(a, n):
    for i in range(n):
        for j in range(n):
            print(a[i][j], end=' ')
        print('')

def run(d, t):
    a_2D = []
    b_2D = []

    for i in range(d):
        row = []
        for j in range(d):
            row.append(random.randint(0, 3))
        a_2D.append(row)

    for i in range(d):
        row = []
        for j in range(d):
            row.append(random.randint(0, 3))
        b_2D.append(row)

    global THRESHOLD
    global BASE_CASE
    THRESHOLD = t

    a_2D, BASE_CASE = pad(a_2D)
    b_2D, BASE_CASE = pad(b_2D)

    # Size of the padded matrix
    padded_d = len(a_2D)

    a_1D = []
    b_1D = []

    # print("\nMATRIX A")
    # print_diagonal(a_2D, d)
    #
    # print("\nMATRIX B")
    # print_diagonal(b_2D, d)

    z_order(a_1D, a_2D, 0, 0, padded_d, m=BASE_CASE)
    z_order(b_1D, b_2D, 0, 0, padded_d, m=BASE_CASE)

    c = strassen(a_1D, b_1D, padded_d)


    # c = [0] * len(a_1D)
    # standard_multiply(a_1D, b_1D, c, n=get_size(a_1D), m=get_size(a_1D), idx=0)

    # c = [0 for i in range(len(a_1D))]
    # standard_multiply(a_1D, b_1D, c, padded_d, m=BASE_CASE, idx=0)
    # c_ordered = [[] for i in range(padded_d)]
    # c_ordered = reorder(c, c_ordered, padded_d, m=BASE_CASE, row=0)

    # c_ordered = [[] for i in range(padded_d)]
    # c_ordered = reorder(c, c_ordered, padded_d, m=BASE_CASE, row=0)

    # print("\nMATRIX C")
    # print_diagonal(c_ordered, d)

    return padded_d


def count_triangles(probabilities):

    for p in probabilities:
        d = 1024

        global THRESHOLD
        global BASE_CASE

        THRESHOLD = 31
        BASE_CASE = 16

        a_2D = [[0 for j in range(d)] for i in range(d)]

        for i in range(d-1):
            for j in range(i+1, d):
                if random.uniform(0, 1) <= p:
                    a_2D[i][j] = 1
                    a_2D[j][i] = 1

        a_1D = []

        z_order(a_1D, a_2D, 0, 0, d, m=BASE_CASE)

        c = strassen(a_1D, strassen(a_1D, a_1D, d), d)

        c_ordered = [[] for i in range(d)]
        c_ordered = reorder(c, c_ordered, d, m=BASE_CASE, row=0)

        sum = 0
        for i in range(d):
            sum += c_ordered[i][i]

        num_of_triangles = sum / 6
        expectation = d*(d-1)*(d-2)*(p**3) / 6
        print(f'{num_of_triangles}, {expectation}')


def run_mp(sub_tuples):
    for single_tuple in sub_tuples:
        t = single_tuple[0]
        d = single_tuple[1]
        start = time.time()
        padded_d = run(d, t)
        end = time.time()
        print(f'{d},{t}, {padded_d}, {BASE_CASE}, {end - start}')


if __name__ == '__main__':
    # input_tuples = [[(64, 1024), (8, 1024), (64, 1024), (8, 1024), (16, 1024), (8, 1024)],
    #                 [(32, 1024), (16, 1024), (32, 1024), (16, 1024), (8, 1024), (16, 1024)],
    #                 [(16, 1024), (32, 1024), (8, 1024), (8, 1024), (16, 1024), (8, 1024)],
    #                 [(8, 1024), (64, 1024), (16, 1024), (16, 1024), (8, 1024), (16, 1024)]]
    # sub_tuples = []

    # visited_bc = set()
    # for i in range(9, 83):
    #
    #     if i in visited_bc:
    #         continue
    #
    #     size = i
    #
    #     while size <= 1024:
    #         size = size * 2
    #
    #     bc = i
    #
    #     while bc <= 167:
    #         sub_tuples.append((bc, size))
    #         visited_bc.add(bc)
    #         bc = bc * 2
    #
    #         if len(sub_tuples) == 30:
    #             input_tuples.append(sub_tuples)
    #             sub_tuples = []

    # # Multiprocessing
    # pool = mp.Pool(processes=len(input_tuples))
    # output_async = [pool.apply_async(func=run_mp, args=(sub_tuples,)) for sub_tuples in input_tuples]
    # output = [x.get() for x in output_async]

    # dim = [1920]
    # thresholds = [15, 30, 15, 30]
    # for d in dim:
    #     for t in thresholds:
    #         start = time.time()
    #         padded_d = run(d, t)
    #         end = time.time()
    #         print(f'{d},{t}, {padded_d}, {BASE_CASE}, {end - start}')


    # for i in range(500, 2001, 100):
    #     start = time.time()
    #     d = i
    #     t = i
    #     run(d, t)
    #     end = time.time()
    #     print(f'{d},{t},{end - start}')

    probs = [0.01, 0.02, 0.03, 0.04, 0.05]
    probabilities = []

    for p in probs:
        sub_list = []
        for i in range(10):
            sub_list.append(p)
        probabilities.append(sub_list)

    # Multiprocessing
    pool = mp.Pool(processes=len(probabilities))
    output_async = [pool.apply_async(func=count_triangles, args=(sub_list,)) for sub_list in probabilities]
    output = [x.get() for x in output_async]

    print('DONE!')


