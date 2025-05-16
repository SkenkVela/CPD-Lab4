#include <mpi.h>
#include <stdio.h>
int main(int argc, char *argv[])
{
    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    printf("Oi do host %d de %d\n", id, p);
    MPI_Finalize();

    return 0;
}


