// largura_banda.c
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MESSAGE_SIZE_MB 100
#define MB (1024 * 1024)

int main(int argc, char *argv[]) {
    int rank;
    char *buffer;
    double start, end, bandwidth;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t size = MESSAGE_SIZE_MB * MB;
    buffer = (char *) malloc(size);
    if (!buffer) {
        if (rank == 0) printf("Falha ao alocar memória.\n");
        MPI_Finalize();
        return 1;
    }

    // Inicializa o buffer com algo (não obrigatório, mas bom)
    for (size_t i = 0; i < size; i++) {
        buffer[i] = 'A';
    }

    if (rank == 0) {
        start = MPI_Wtime();
        MPI_Send(buffer, size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(buffer, 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status);
        end = MPI_Wtime();

        double time_sec = end - start;
        bandwidth = (size / (1024.0 * 1024.0)) / time_sec; // MB/s

        printf("Largura de banda: %.2f MB/s, Tempo: %.6f segundos\n", bandwidth, time_sec);
    } else if (rank == 1) {
        MPI_Recv(buffer, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Send(buffer, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    free(buffer);
    MPI_Finalize();
    return 0;
}
 