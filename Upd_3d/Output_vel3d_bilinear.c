/* 可替换原 Output_vel3d；需要项目原有的 count_file、Arraydouble2、
 * rdfile、Freedouble2 声明。剖面输入为 x z Vp Vs rho 的完整规则网格。
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int compare_double(const void *left, const void *right)
{
    double a = *(const double *)left, b = *(const double *)right;
    return (a > b) - (a < b);
}

static int unique_sorted(double *axis, int length)
{
    int i, count = 0;
    qsort(axis, (size_t)length, sizeof(*axis), compare_double);
    for (i = 0; i < length; i++) {
        if (!count || axis[i] != axis[count - 1]) axis[count++] = axis[i];
    }
    return count;
}

/* 在有序轴中返回左端网格编号和单元内比例；超出范围时贴到边界。 */
static int cell_position(const double *axis, int count, double value, double *fraction)
{
    int low = 0, high = count - 1, middle;
    if (value <= axis[0]) { *fraction = 0.0; return 0; }
    if (value >= axis[high]) { *fraction = 1.0; return high - 1; }
    while (high - low > 1) {
        middle = low + (high - low) / 2;
        if (axis[middle] <= value) low = middle;
        else high = middle;
    }
    *fraction = (value - axis[low]) / (axis[low + 1] - axis[low]);
    return low;
}

void Output_vel3d(short num, long node, long elem, long **Elem_array, int Ns, int Nr,
                  float *rho, float *lamda, float *mu, double **node_array)
{
    int n, numVel, nx, nz, ix, iz, component;
    long i, j, node_id, corner[4];
    size_t k, total;
    double sx, sz, tx, tz, x, z, result[3];
    double **arrayVel = NULL, *xaxis = NULL, *zaxis = NULL, *grid = NULL;
    char fileVelIt[256], fileVelOt[256];
    FILE *fp = NULL;

    (void)rho;
    (void)lamda;
    (void)mu;

    if (num < 2 || node <= 0 || elem <= 0 || !Elem_array || !node_array) {
        fprintf(stderr, "Output_vel3d: 网格参数无效\n");
        exit(EXIT_FAILURE);
    }
    n = snprintf(fileVelIt, sizeof fileVelIt,
                 "../Output/Model/profile2d/profile2d_s%d_r%d.txt", Ns, Nr);
    if (n < 0 || (size_t)n >= sizeof fileVelIt) {
        fprintf(stderr, "Output_vel3d: 输入路径过长\n");
        exit(EXIT_FAILURE);
    }
    n = snprintf(fileVelOt, sizeof fileVelOt,
                 "../Output/Model/Forver/vel2d_s%d_r%d.txt", Ns, Nr);
    if (n < 0 || (size_t)n >= sizeof fileVelOt) {
        fprintf(stderr, "Output_vel3d: 输出路径过长\n");
        exit(EXIT_FAILURE);
    }

    numVel = count_file(fileVelIt);
    if (numVel <= 0) {
        fprintf(stderr, "Output_vel3d: 无法读取或剖面为空：%s\n", fileVelIt);
        exit(EXIT_FAILURE);
    }
    arrayVel = Arraydouble2(numVel, 5);
    xaxis = malloc((size_t)numVel * sizeof *xaxis);
    zaxis = malloc((size_t)numVel * sizeof *zaxis);
    if (!arrayVel || !xaxis || !zaxis) {
        fprintf(stderr, "Output_vel3d: 内存分配失败\n");
        exit(EXIT_FAILURE);
    }

    rdfile(fileVelIt, numVel, 5, arrayVel);

    for (i = 0; i < numVel; i++) {
        for (j = 0; j < 5; j++) {
            if (!isfinite(arrayVel[i][j])) {
                fprintf(stderr, "Output_vel3d: 剖面含无效数值（第 %ld 行）\n", i + 1);
                exit(EXIT_FAILURE);
            }
        }
        xaxis[i] = arrayVel[i][0];
        zaxis[i] = arrayVel[i][1];
    }
    
    nx = unique_sorted(xaxis, numVel);
    nz = unique_sorted(zaxis, numVel);
    total = (size_t)nx * (size_t)nz;
    if (nx < 2 || nz < 2 || total != (size_t)numVel || total > SIZE_MAX / (3 * sizeof *grid)) {
        fprintf(stderr, "Output_vel3d: 剖面必须是至少 2x2 的完整 x-z 规则网格\n");  exit(EXIT_FAILURE);
    }
    grid = malloc(3 * total * sizeof *grid);
    if (!grid) {
        fprintf(stderr, "Output_vel3d: 网格分配失败\n");  exit(EXIT_FAILURE);
    }
    for (k = 0; k < 3 * total; k++) grid[k] = NAN;
    /* 不依赖输入文件的行排列，转为 x 外层、z 内层。 */
    for (i = 0; i < numVel; i++) {
        const double *px = bsearch(&arrayVel[i][0], xaxis, (size_t)nx, sizeof *xaxis, compare_double);
        const double *pz = bsearch(&arrayVel[i][1], zaxis, (size_t)nz, sizeof *zaxis, compare_double);
        if (!px || !pz) {
            fprintf(stderr, "Output_vel3d: 剖面坐标不在网格轴上\n");  exit(EXIT_FAILURE);
        }
        k = ((size_t)(px - xaxis) * (size_t)nz + (size_t)(pz - zaxis)) * 3;
        if (isfinite(grid[k])) {
            fprintf(stderr, "Output_vel3d: 剖面中有重复的 (x,z) 点\n");  exit(EXIT_FAILURE);
        }
        for (component = 0; component < 3; component++)
            grid[k + (size_t)component] = arrayVel[i][component + 2];
    }

    for (k = 0; k < 3 * total; k++) {
        if (!isfinite(grid[k])) {
            fprintf(stderr, "Output_vel3d: 剖面网格缺少采样点\n");  exit(EXIT_FAILURE);
        }
    }

    fp = fopen(fileVelOt, "w");
    if (!fp) { perror(fileVelOt);  exit(EXIT_FAILURE); }
    corner[0] = 0;
    corner[1] = (long)num * (num - 1);
    corner[2] = (long)num * num - 1;
    corner[3] = num - 1;

    for (i = 0; i < elem; i++) {
        sx = sz = 0.0;
        for (j = 0; j < 4; j++) {
            node_id = Elem_array[i][corner[j]] - 1;
            if (node_id < 0 || node_id >= node) {
                fprintf(stderr, "Output_vel3d: 第 %ld 个单元的节点编号无效\n", i + 1);  exit(EXIT_FAILURE);
            }
            sx += node_array[node_id][0];
            sz += node_array[node_id][1];
        }
        x = sx / 4.0;
        z = sz / 4.0;
        ix = cell_position(xaxis, nx, x, &tx);
        iz = cell_position(zaxis, nz, z, &tz);
        for (component = 0; component < 3; component++) {
            size_t a = ((size_t)ix * (size_t)nz + (size_t)iz) * 3 + (size_t)component;
            size_t b = a + (size_t)nz * 3;
            result[component] = (1 - tx) * (1 - tz) * grid[a]
                              + tx * (1 - tz) * grid[b]
                              + (1 - tx) * tz * grid[a + 3]
                              + tx * tz * grid[b + 3];
        }
        if (fprintf(fp, "%.4f %.4f %.4f %.4f %.4f\n",
                    x, z, result[0], result[1], result[2]) < 0) {
            perror(fileVelOt);
            exit(EXIT_FAILURE);
        }
    }
    if (fclose(fp) != 0) { perror(fileVelOt); exit(EXIT_FAILURE); }
    Freedouble2(numVel, arrayVel);
    free(xaxis); free(zaxis); free(grid);
}
