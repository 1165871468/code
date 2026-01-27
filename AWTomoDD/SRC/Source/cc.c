#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // 每个进程的本地二维数组（3x4）
    int rows = 3, cols = 4;
    double local_data[rows][cols];
    
    // 初始化本地数据
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            local_data[i][j] = rank * 100.0 + i * cols + j + 0.5;
        }
    }
    
    // 根进程准备接收缓冲区
    double* gathered_data = NULL;
    int* recvcounts = NULL;
    int* displs = NULL;
    
    if (rank == 0) {
        gathered_data = (double*)malloc(size * rows * cols * sizeof(double));
        recvcounts = (int*)malloc(size * sizeof(int));
        displs = (int*)malloc(size * sizeof(int));

        for (int i = 0; i < size; i++) {
            recvcounts[i] = rows * cols;  // 每个进程发送的数据量
            displs[i] = i * rows * cols;  // 在接收缓冲区中的偏移量
        }
    }

    // 收集所有进程的数据
    MPI_Gatherv(local_data,           // 发送缓冲区
                rows * cols,          // 发送数据量
                MPI_DOUBLE,           // 数据类型改为 MPI_DOUBLE
                gathered_data,        // 接收缓冲区
                recvcounts,           // 每个进程接收的数据量
                displs,               // 接收偏移量
                MPI_DOUBLE,           // 数据类型改为 MPI_DOUBLE
                0,                    // 根进程
                MPI_COMM_WORLD);
    
    // 根进程输出结果
    if (rank == 0) {
        printf("Gathered 2D arrays (double):\n");
        for (int p = 0; p < size; p++) {
            printf("Process %d:\n", p);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    int idx = p * rows * cols + i * cols + j;
                    printf("%8.2f ", gathered_data[idx]);
                }
                printf("\n");
            }
            printf("\n");
        }
        free(gathered_data);
        free(recvcounts);
        free(displs);
    }
    
    MPI_Finalize();
    return 0;
}

viod MPI_GetData2d(int size, int rows, int cols, double **arrayIn, double **arrayOut)
{
    int *recvcounts, *displs;
    double *gathered_data;

    gathered_data = (double*)malloc(size * rows * cols * sizeof(double));
    recvcounts = (int*)malloc(size * sizeof(int));
    displs = (int*)malloc(size * sizeof(int));


    MPI_Gatherv(arrayIn,           // 发送缓冲区
                rows * cols,          // 发送数据量
                MPI_DOUBLE,           // 数据类型改为 MPI_DOUBLE
                gathered_data,        // 接收缓冲区
                recvcounts,           // 每个进程接收的数据量
                displs,               // 接收偏移量
                MPI_DOUBLE,           // 数据类型改为 MPI_DOUBLE
                0,                    // 根进程
                MPI_COMM_WORLD);

    free(gathered_data);
    free(recvcounts);
    free(displs);
}