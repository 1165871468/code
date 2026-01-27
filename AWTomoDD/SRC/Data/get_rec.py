import math
import numpy as np
import time
import matplotlib.pyplot as plt

def find_nest(array, vv):
    id = np.abs(array - vv).argmin()

    return id

def getnum(flnm):
    id1 = 0

    with open(flnm) as file:
        for line in file:
            line = line.replace('\n', '')
            cc = line.split()

            if id1 == 1:
                num_ele = int(cc[0])
                num_nod = int(cc[1])
                id1 = 2
                break

            if cc[0] == 'number':
                id1 = 1

    print('number of elements: ', num_ele)
    print('number of nodes: ', num_nod)

    return num_nod, num_ele

def getnods(flnm, num_nod):
    data1 = np.zeros((num_nod, 3))

    id1 = 0; id2 = 0

    with open(flnm) as file:
        for line in file:
            line = line.replace('\n', '')
            cc = line.split()

            if id1 == 1 and cc[0] != 'End':
                data1[id2][0] = float(cc[0])
                data1[id2][1] = float(cc[1])
                data1[id2][2] = float(cc[2])

                id2 = id2 + 1
            elif cc[0] == 'End':
                break

            if cc[0] == 'coordinates':
                id1 = 1

    return data1

def geteles(flnm, num_ele):
    data2 = np.zeros((num_ele, 6))

    id1 = 0; id2 = 0

    with open(flnm) as file:
        for line in file:
            line = line.replace('\n', '')
            cc = line.split()

            if id1 == 1 and cc[0] != 'End':
                data2[id2][0] = float(cc[0])
                data2[id2][1] = float(cc[1])
                data2[id2][2] = float(cc[2])
                data2[id2][3] = float(cc[3])
                data2[id2][4] = float(cc[4])
                data2[id2][5] = float(cc[5])

                id2 = id2 + 1
            elif id1 == 1 and cc[0] == 'End':
                break

            if cc[0] == 'elements':
                id1 = 1

    return data2

def find_neast_point(points, x0, y0):
    dis = np.sqrt( (points[:, 0] - x0)**2 + (points[:, 1] - y0)**2 )

    id = np.argmin(dis)

    point = points[id]

    return point


if __name__ == '__main__':

    num_rec = 21

    filnm1 = 'topo3.txt'
    filnm2 = 'receiver'
    filnm3 = 'Mesh_test.txt'

    am = 4000

    data = np.loadtxt(filnm1)

    x = data[:, 0]
    y = data[:, 1]

    count1 = len(x)
    num_rec = num_rec - 1

    dx1 = x[-1] / num_rec
    dx2 = math.floor(count1 / num_rec)

    rec = np.zeros((num_rec + 1, 2))

    for i in range(num_rec + 1):
        xx = i * dx1

        id = find_nest(x, xx)

        zz = math.floor(y[id])

        rec[i][0] = xx
        rec[i][1] = zz + 200


    num_nod, num_ele = getnum(filnm3)

    data1 = getnods(filnm3, num_nod)

    xx = data1[:, 1]
    yy = -data1[:, 2] + am

    points = np.column_stack((xx, yy))

    for i in range(num_rec + 1):
        x0 = rec[i][0]
        y0 = rec[i][1]

        rec[i] = find_neast_point(points, x0, y0)

    np.savetxt(filnm2, rec, fmt = '%.1f %10.0f')
