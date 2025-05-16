#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ARRAY_SIZE_MB 100
#define MB (1024 * 1024)

int main(int argc, char *argv[]) {
    int rank, size;
    size_t array_size = ARRAY_SIZE_MB * MB;
    char *array = NULL;
    double start, end, time_bcast, time_send;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Aloca o array só no root (0)
    if (rank == 0) {
        array = (char *) malloc(array_size);
        if (!array) {
            printf("Erro na alocação de memória\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Inicializa o array
        for (size_t i = 0; i < array_size; i++) {
            array[i] = 'A';
        }
    } else {
        // Nos outros processos, aloca para receber
        array = (char *) malloc(array_size);
        if (!array) {
            printf("Erro na alocação de memória\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // --- MPI_Bcast ---
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    MPI_Bcast(array, array_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    time_bcast = end - start;

    // --- MPI_Send individual ---
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    if (rank == 0) {
        for (int dest = 1; dest < size; dest++) {
            MPI_Send(array, array_size, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(array, array_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    time_send = end - start;

    if (rank == 0) {
        printf("Tempo MPI_Bcast: %.6f segundos\n", time_bcast);
        printf("Tempo MPI_Send individual: %.6f segundos\n", time_send);
        printf("Número de processos: %d\n", size);
        printf("Tamanho do array: %d MB\n", ARRAY_SIZE_MB);
    }

    free(array);
    MPI_Finalize();
    return 0;
}
