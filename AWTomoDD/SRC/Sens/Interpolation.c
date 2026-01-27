#include "Common.h"
#include "Interpolation.h"

int Get_num_xz(int count, double x1, double x2, double xx)
{
    int i, id;
    double dx, xx1, *array_x;

    dx = (x2 - x1)/(count - 1);

    id = -1;
    for (i = 0; i < count; i++)
    {
        xx1 = i*dx + x1;
        if (xx >= xx1 && xx < (xx1+dx))
        {
            id = i;
            break;
        }
    }

    if (id == count && xx > ((count - 1)*dx + x1))
    {
        id = -1;
    }

    if (id == -1)
    {
        printf("Can not find margin %.1f -- %.1f Vs %.1f\n", xx, x1, x2);
    }
    
    return id;
}

int Nearest_point1d(double target, double *points, int count)
{
    int i, id;
    double flgd1, distance, min_dis = INFINITY;

    for (i = 0; i < count; i++)
    {
        distance = fabs(target - points[i]);

        if (distance < min_dis)
        {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point2d(double *target, double **points, int count)
{
    int i, id;
    double flgd1, distance, min_dis = INFINITY;

    id = -12345;
    for (i = 0; i < count; i++)
    {
        distance = sqrt(pow(target[0] - points[i][0], 2) + pow(target[1] - points[i][1], 2));

        if (distance < min_dis)
        {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point3d(double *target, double **points, int count)
{
    int id;
    double flgd1, min_dis = INFINITY;

    for (int i = 0; i < count; i++)
    {
        double distance = sqrt(pow(target[0]-points[i][0], 2) + pow(target[1]-points[i][1], 2) + 1.0*pow(target[2]-points[i][2], 2));

        if (distance < min_dis)
        {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point2d_lu(double *target, double **points, int count)
{
    int id;
    double flgd1, min_dis = INFINITY;

    id = -12345;
    for (int i = 0; i < count; i++) {
        double distance = sqrt(pow(target[0] - points[i][0], 2) + pow(target[1] - points[i][1], 2));

        if (distance < min_dis && points[i][0] <= target[0] && points[i][1] >= target[1]) {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point2d_ru(double *target, double **points, int count)
{
    int id;
    double flgd1, min_dis = INFINITY;

    id = -12345;
    for (int i = 0; i < count; i++) {
        double distance = sqrt(pow(target[0]-points[i][0], 2) + pow(target[1]-points[i][1], 2));

        if (distance < min_dis && points[i][0] >= target[0] && points[i][1] >= target[1]) {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point2d_ld(double *target, double **points, int count)
{
    int id;
    double flgd1, min_dis = INFINITY;

    id = -12345;
    for (int i = 0; i < count; i++) {
        double distance = sqrt(pow(target[0]-points[i][0], 2) + pow(target[1]-points[i][1], 2));

        if (distance < min_dis && points[i][0] <= target[0] && points[i][1] <= target[1]) {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

int Nearest_point2d_rd(double *target, double **points, int count)
{
    int id;
    double flgd1, min_dis = INFINITY;

    id = -12345;
    for (int i = 0; i < count; i++) {
        double distance = sqrt(pow(target[0]-points[i][0], 2) + pow(target[1]-points[i][1], 2));

        if (distance < min_dis && points[i][0] >= target[0] && points[i][1] <= target[1]) {
            min_dis = distance;
            id = i;
        }
    }
    return id;
}

// 计算二维向量的交叉乘积
double Cross(Vec2 a, Vec2 b)
{
    return a.x*b.y - a.y*b.x;
}


double Interplot_1d(double xx[3], double yy[3])
{
    int i;
    double aa, bb, coff, ff;

    ff = 0.0;
    if ((xx[0] - xx[1]) == 0.0)
    {
        xx[1] = xx[0] + 0.1;
    }

    for (i = 0; i < 2; i++)
    {
        aa = xx[2] - xx[i];
        bb = xx[1] - xx[0];
        coff = 1.0 - fabs(aa/bb);
        ff += coff*yy[i];
    }
    return ff;
}

double Interplot_2d(double xx[3], double yy[3], double zz[2][2])
{
    int i, j;
    double aa, bb, coff1, coff2, ff;

    ff = 0.0;

    if ((xx[1] - xx[0]) == 0.0)
    {
        xx[1] = xx[0] + 0.1;
    }
    if ((yy[1] - yy[0]) == 0.0)
    {
        yy[1] = yy[0] + 0.1;
    }

    for (i = 0; i < 2; i++)
    {
        aa = xx[2] - xx[i];
        bb = xx[1] - xx[0];
        coff1 = 1.0 - fabs(aa/bb);

        for (j = 0; j < 2; j++)
        {
            aa = yy[2] - yy[j];
            bb = yy[1] - yy[0];
            coff2 = 1.0 - fabs(aa/bb);
            ff += coff1*coff2*zz[i][j];
        }
    }
    return ff;
}

double Interplot_2din(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d, double fa, double fb, double fc, double fd)
{
    double u, v, result;
    Vec2 uv = Invbilinear(p, a, b, c, d);

    if (uv.x < 0.0 || uv.y < 0.0) {
        return -1.0;
    }
    u = uv.x;
    v = uv.y;
    result = (1.0 - u)*(1.0 -v)*fa + u*(1.0 -v)*fb + u*v*fc + (1.0 - u)*v*fd;

    return result;
}

double Interplot_3d(double xx[3], double yy[3], double zz[3], double vv[2][2][2])
{
    int i, j, k;
    double aa, bb, coff1, coff2, coff3, ff;

    ff = 0.0;

    if ((xx[1] - xx[0]) == 0.0)
    {
        xx[1] = xx[0] + 0.1;
    }
    if ((yy[1] - yy[0]) == 0.0)
    {
        yy[1] = yy[0] + 0.1;
    }
    if ((zz[1] - zz[0]) == 0.0)
    {
        zz[1] = zz[0] + 0.1;
    }

    for (i = 0; i < 2; i++)
    {
        aa = xx[2] - xx[i];
        bb = xx[1] - xx[0];
        coff1 = 1.0 - fabs(aa/bb);

        for (j = 0; j < 2; j++)
        {
            aa = yy[2] - yy[j];
            bb = yy[1] - yy[0];
            coff2 = 1.0 - fabs(aa/bb);

            for (k = 0; k < 2; k++)
            {
                aa = zz[2] - zz[k];
                bb = zz[1] - zz[0];
                coff3 = 1.0 - fabs(aa/bb);
                ff += coff1*coff2*coff3*vv[i][j][k];
            }
        }
    }

    return ff;
}

Vec2 Invbilinear(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d)
{
    Vec2 res = {-1.0, -1.0};
    Vec2 e = {b.x - a.x, b.y - a.y};
    Vec2 f = {d.x - a.x, d.y - a.y};
    Vec2 g = {a.x - b.x + c.x - d.x, a.y - b.y + c.y - d.y};
    Vec2 h = {p.x - a.x, p.y - a.y};

    double k2 = Cross(g, f);
    double k1 = Cross(e, f) + Cross(h, g);
    double k0 = Cross(h, e);

    if (fabs(k2) < 0.001f) {
        res.x = (h.x*k1 + f.x*k0)/(e.x*k1 - g.x*k0);
        res.y = -k0/k1;
    }
    else {
        double w = k1*k1 - 4.0*k0*k2;

        if (w < 0.0) {
            return res;
        }
        w = sqrt(w);

        double ik2 = 0.5 / k2;
        double v = (-k1 - w)*ik2;
        double u = (h.x - f.x*v) / (e.x + g.x*v);

        if (u < 0.0 || u > 1.0 || v < 0.0 || v > 1.0) {
            v = (-k1 + w)*ik2;
            u = (h.x - f.x*v)/(e.x + g.x*v);
        }

        res.x = u;
        res.y = v;
    }

    return res;
}

void Find_node(int numData, double **arrayData, double *tarPoints, double arrayPoints[][3])
{
    int id1, id2, id3, id4;
    double coff, ff, target[2];

    target[0] = tarPoints[0];
    target[1] = tarPoints[1];

    id1 = Nearest_point2d_lu(target, arrayData, numData);
    id2 = Nearest_point2d_ru(target, arrayData, numData);
    id3 = Nearest_point2d_ld(target, arrayData, numData);
    id4 = Nearest_point2d_rd(target, arrayData, numData);

    arrayPoints[0][0] = arrayData[id1][0];  arrayPoints[0][1] = arrayData[id1][1];  arrayPoints[0][2] = arrayData[id1][2];
    arrayPoints[1][0] = arrayData[id2][0];  arrayPoints[1][1] = arrayData[id2][1];  arrayPoints[1][2] = arrayData[id2][2];
    arrayPoints[2][0] = arrayData[id3][0];  arrayPoints[2][1] = arrayData[id3][1];  arrayPoints[2][2] = arrayData[id3][2];
    arrayPoints[3][0] = arrayData[id4][0];  arrayPoints[3][1] = arrayData[id4][1];  arrayPoints[3][2] = arrayData[id4][2];
}

void Fine_node_cube(int numY, int numZ, int idx1, int idy1, int idz1, int idx2, int idy2, int idz2, 
    double **arrayV, double xx[3], double yy[3], double zz[3], double vv[2][2][2])
{
    int i, j, id[8];

    id[0] = idz1 + numZ*idy1 + numY*numZ*idx1;
    id[1] = idz1 + numZ*idy1 + numY*numZ*idx2;
    id[2] = idz2 + numZ*idy1 + numY*numZ*idx1;
    id[3] = idz2 + numZ*idy1 + numY*numZ*idx2;
    id[4] = idz1 + numZ*idy2 + numY*numZ*idx1;
    id[5] = idz1 + numZ*idy2 + numY*numZ*idx2;
    id[6] = idz2 + numZ*idy2 + numY*numZ*idx1;
    id[7] = idz2 + numZ*idy2 + numY*numZ*idx2;

    xx[0] = arrayV[id[0]][0];         yy[0] = arrayV[id[0]][1];         zz[0] = arrayV[id[0]][2];
    xx[1] = arrayV[id[7]][0];         yy[1] = arrayV[id[7]][1];         zz[1] = arrayV[id[7]][2];

    vv[0][0][0] = arrayV[id[0]][3];   vv[0][1][0] = arrayV[id[4]][3];
    vv[1][0][0] = arrayV[id[1]][3];   vv[1][1][0] = arrayV[id[5]][3];
    vv[0][0][1] = arrayV[id[2]][3];   vv[0][1][1] = arrayV[id[6]][3];
    vv[1][0][1] = arrayV[id[3]][3];   vv[1][1][1] = arrayV[id[7]][3];
}

int Fine_node_square(int numData, int id[4], double **arrayData, double *locNode)
{
    int i, j, flgi;

    id[0] = Nearest_point2d_lu(locNode, arrayData, numData);
    id[1] = Nearest_point2d_ru(locNode, arrayData, numData);
    id[2] = Nearest_point2d_rd(locNode, arrayData, numData);
    id[3] = Nearest_point2d_ld(locNode, arrayData, numData);

    flgi = 1;
    for (i = 0; i < 4; i++) {
        if (id[i] == -12345) {
            flgi = -1;
        }
    }

    return flgi;
}
