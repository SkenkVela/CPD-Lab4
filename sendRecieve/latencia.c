// latencia.c
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, p;
    int msg = 0;
    double start, end, latency;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (p < 2) {
        if (rank == 0) printf("Este teste requer pelo menos 2 processos.\n");
        MPI_Finalize();
        return 1;
    }

    // Apenas o processo 0 e 1 participam
    if (rank == 0) {
        start = MPI_Wtime(); // Início da medição
        MPI_Send(&msg, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&msg, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        end = MPI_Wtime(); // Fim da medição

        latency = (end - start) / 2.0;
        printf("Latência (ida ou volta) = %.9f segundos = %.3f microssegundos\n", latency, latency * 1e6);
    }
    else if (rank == 1) {
        MPI_Recv(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
