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


def multiply(a, b, c, n, k):

    if n == k:
        for i in range(k):
            r = 0
            for j in range(k):
                r += a[i] * b[i*k]

    else:
        n = n // 2

        p = n // 4

        a00 = a[:p]
        a01 = a[p:n]
        a10 = a[n:n+p]
        a11 = a[n+p:]

        b00 = b[:p]
        b01 = b[p:n]
        b10 = b[n:n + p]
        b11 = b[n + p:]

        c00 = c[:p]
        c01 = c[p:n]
        c10 = c[n:n + p]
        c11 = c[n + p:]

        multiply(a00, b00, c00, n, k)
        multiply(a01, b10, c00, n, k)
        multiply(a00, b01, c01, n, k)
        multiply(a01, b11, c01, n, k)





if __name__ == '__main__':
    a = [[0, 1, 2, 2, 1, 2],
         [2, 1, 0, 1, 1, 2],
         [2, 2, 2, 1, 0, 2],
         [1, 0, 2, 2, 0, 1],
         [2, 1, 0, 1, 0, 1],
         [2, 2, 1, 1, 0, 1]]
    a_morton = []

    b = [[0, 1, 2, 1, 0, 0],
         [1, 0, 2, 0, 1, 2],
         [1, 0, 0, 1, 1, 2],
         [1, 2, 1, 2, 2, 2],
         [0, 2, 1, 2, 0, 0],
         [0, 1, 1, 2, 2, 2]]
    b_morton = []

    z_order(a_morton, a, 0, 0, len(a), k=3)
    z_order(b_morton, b, 0, 0, len(b), k=3)

    c = [0 for i in range(len(a_morton))]
    c = multiply(a_morton, b_morton, c)

    d = len(a)
    for i in range(d):
        for j in range(d):
            print(c[i][j], end=' ')
        print('')

