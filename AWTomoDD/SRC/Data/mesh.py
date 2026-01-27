import math
import numpy as np
import time

def Test_mesh(flname, m, n):
    m = m + 1
    n = n + 1
    det = 1000.0

    x = np.zeros(m)
    y = np.zeros(n)

    for i in range(m):
        x[i] = i * det

    for i in range(n):
        y[i] = i * det

    fp = open(flname, 'w')

    for i in range(m):
        for j in range(n):
            print(f'{x[i]:.1f} {y[j]:.1f}', file = fp); fp.flush

    fp.close()

def Get_elem_node(flname1, flname2, m, n, x0, x1):
    num_ele = n * m

    m = m + 1
    n = n + 1

    num_node = m * n

    data = np.loadtxt(flname1)

    count = len(data)

    if num_node != count:
        print("The data is not marching!!")
        exit(1)

    ele = np.zeros((num_ele + 1, 4))
    node = np.zeros((num_node + 1, 3))

    fp1 = open(flname2, 'w')

    print('number of elements & nodes & ndime:', file = fp1); fp1.flush()
    print(f'{num_ele:d} {num_node:d} {2:d}', file = fp1); fp1.flush()
    print('coordinates', file = fp1); fp1.flush()

    for i in range(num_node):
        id = i + 1

        Col = id % m

        if Col == 0:
            Col = m
            Row = id // m - 1
        else:
            Row = id // m

        num = int(Row + (Col - 1) * n)

        node[i][0] = i + 1
        node[i][1] = data[num][0]
        node[i][2] = data[num][1]

        print(f'{node[i][0]:.0f} {node[i][1]:.1f} {node[i][2]:.1f}', file = fp1); fp1.flush()

    print('End coordinates', file = fp1); fp1.flush()
    print('elements', file = fp1); fp1.flush()

    M = m
    m = m - 1

    for i in range(num_ele):
        id = i + 1

        col_m = id % m

        if col_m == 0:
            ln_m = id // m - 1

            col_m = m
        else:
            ln_m = id // m

        tmp1 = ln_m * M + col_m
        tmp2 = tmp1 + 1
        tmp3 = tmp2 + M
        tmp4 = tmp1 + M

        # ele[i][0] = tmp2
        # ele[i][1] = tmp1
        # ele[i][2] = tmp4
        # ele[i][3] = tmp3

        ele[i][0] = tmp3
        ele[i][1] = tmp4
        ele[i][2] = tmp1
        ele[i][3] = tmp2

        xx = (node[tmp1-1][1] + node[tmp2-1][1] + node[tmp3-1][1] + node[tmp4-1][1]) / 4

        if xx > x0 and xx < x1:
            # print(xx, x0, x1)
            # print(f'{id:.0f} {ele[i][0]:.0f} {ele[i][1]:.0f} {ele[i][2]:.0f} {ele[i][3]:.0f} {1:d}')
            # print('+++++++++++++++++++++++++++++++=')
            print(f'{id:.0f} {ele[i][0]:.0f} {ele[i][1]:.0f} {ele[i][2]:.0f} {ele[i][3]:.0f} {1:d}', file = fp1); fp1.flush()


    for i in range(num_ele):
        id = i + 1

        col_m = id % m

        if col_m == 0:
            ln_m = id // m - 1

            col_m = m
        else:
            ln_m = id // m

        tmp1 = ln_m * M + col_m
        tmp2 = tmp1 + 1
        tmp3 = tmp2 + M
        tmp4 = tmp1 + M

        ele[i][0] = tmp3
        ele[i][1] = tmp4
        ele[i][2] = tmp1
        ele[i][3] = tmp2

        xx = (node[tmp1-1][1] + node[tmp2-1][1] + node[tmp3-1][1] + node[tmp4-1][1]) / 4

        if xx < 0 or xx > x1:
            print(f'{id:.0f} {ele[i][0]:.0f} {ele[i][1]:.0f} {ele[i][2]:.0f} {ele[i][3]:.0f} {2:d}', file = fp1); fp1.flush()

    print('End elements', file = fp1); fp1.flush()


def Get_mesh(flnm, m, n, dx, dy):

    for i in range(m):
        for j in range(n):
            continue



if __name__ == '__main__':

    m = 20;    n = 30;     dp = 10
    x0 = 0;    x1 = 40000

    m = m + 2 * dp

    # flname1 = 'Node_test.txt'
    flname1 = 'topo_resple1.txt'
    flname2 = 'Mesh_test.txt'


    # Test_mesh(flname1, m, n)
    Get_elem_node(flname1, flname2, m, n, x0, x1)
