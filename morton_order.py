import math


def z_order(a, b, x, y, n, k):

    if n == k:
        for i in range(k):
            for j in range(k):
                a.append(b[y+i][x+j])
    else:
        n = n // 2

        z_order(a, b, x, y, n, k)
        z_order(a, b, x + n, y, n, k)
        z_order(a, b, x, y + n, n, k)
        z_order(a, b, x + n, y + n, n, k)


def add(a, b):
    c = []
    for i in range(len(a)):
        c.append(a[i] + b[i])

    return c


def subtract(a, b):
    c = []
    for i in range(len(a)):
        c.append(a[i] - b[i])

    return c


def multiply(a, b, c, n, m, idx):

    if n == m:
        for i in range(m):
            for j in range(m):
                r = 0
                for k in range(m):
                    r += a[i * m + k] * b[j + k * m]
                c[idx] += r
                idx += 1

    else:
        p = len(a) // 4
        n = n // 2

        a00 = a[:p]
        a01 = a[p:2*p]
        a10 = a[2*p:3*p]
        a11 = a[3*p:]

        b00 = b[:p]
        b01 = b[p:2*p]
        b10 = b[2*p:3*p]
        b11 = b[3*p:]

        c = multiply(a00, b00, c, n, m, idx)
        c = multiply(a01, b10, c, n, m, idx)
        c = multiply(a00, b01, c, n, m, idx + p)
        c = multiply(a01, b11, c, n, m, idx + p)
        c = multiply(a10, b00, c, n, m, idx + 2*p)
        c = multiply(a11, b10, c, n, m, idx + 2*p)
        c = multiply(a10, b01, c, n, m, idx + 3*p)
        c = multiply(a11, b11, c, n, m, idx + 3*p)

    return c


def reorder(c, d, n, m, row):

    if n == m:
        for i in range(m):
            for j in range(m):
                d[row + i].append(c[i * m + j])

    else:
        p = len(c) // 4
        n = n // 2

        c00 = c[:p]
        c01 = c[p:2*p]
        c10 = c[2*p:3*p]
        c11 = c[3*p:]

        d = reorder(c00, d, n, m, row)

        d = reorder(c01, d, n, m, row)

        d = reorder(c10, d, n, m, row + n)

        d = reorder(c11, d, n, m, row + n)

    return d

def print_matrix(a):
    n = len(a)
    for i in range(n):
        for j in range(n):
            print(d[i][j], end=' ')
        print('')


if __name__ == '__main__':
    a = [[0, 1, 3, 2, 1, 2, 3, 1],
         [2, 1, 0, 3, 1, 2, 3, 1],
         [2, 3, 2, 1, 0, 2, 3, 1],
         [1, 0, 2, 2, 3, 1, 3, 1],
         [3, 1, 0, 1, 0, 1, 3, 1],
         [2, 2, 1, 1, 0, 3, 3, 1],
         [2, 1, 0, 3, 1, 2, 3, 1],
         [2, 3, 2, 1, 0, 2, 3, 1]]
    a_morton = []

    b = [[3, 1, 2, 1, 0, 3, 2, 0],
         [1, 0, 3, 3, 1, 2, 2, 0],
         [1, 3, 0, 1, 1, 3, 2, 0],
         [3, 2, 1, 2, 3, 2, 2, 0],
         [0, 3, 1, 2, 3, 0, 2, 0],
         [3, 1, 1, 3, 2, 2, 2, 0],
         [3, 1, 2, 1, 0, 3, 2, 0],
         [1, 0, 3, 3, 1, 2, 2, 0]]
    b_morton = []

    z_order(a_morton, a, 0, 0, len(a), k=2)
    z_order(b_morton, b, 0, 0, len(b), k=2)

    c = [0 for i in range(len(a_morton))]
    c = multiply(a_morton, b_morton, c, len(a), m=2, idx=0)

    n = int(math.sqrt(len(c)))

    d = [[] for i in range(n)]
    d = reorder(c, d, n, m=2, row=0)
    print_matrix(d)



