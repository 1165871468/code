#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // 补充string.h头文件，用于strlen函数
#include <mpi.h>

int main(int argc, char* argv[]) {
    int rank, size;
    MPI_File file;
    char* filename = "mpi_output.txt";
    char buffer[100];  // 用于存储每个进程要写入的数据

    // 初始化MPI环境
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);




    // 每个进程准备要写入的数据
    // sprintf(buffer, "Hello from process %d of %d\n", rank, size);
    int myLen = snprintf(buffer, sizeof(buffer), "Hello from process %d of %d\n", rank, size);

    MPI_Offset offset = 0;
    MPI_Exscan(&myLen, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // 设置文件指针到正确的位置
    // MPI_File_seek(file, offset, MPI_SEEK_SET);


    // 打开文件，使用MPI_MODE_WRONLY | MPI_MODE_CREATE模式
    // MPI_info参数设为MPI_INFO_NULL表示使用默认设置
    MPI_File_open(MPI_COMM_WORLD, filename, 
                 MPI_MODE_WRONLY | MPI_MODE_CREATE, 
                 MPI_INFO_NULL, &file);
                 
    // 写入数据
    MPI_File_write_at(file, offset, buffer, myLen, MPI_CHAR, MPI_STATUS_IGNORE);

    // 关闭文件
    MPI_File_close(&file);

    // 结束MPI环境
    MPI_Finalize();

    return 0;
}

