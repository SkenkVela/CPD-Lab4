#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

// Implementações simplificadas

int myMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank == root) {
        for (int i = 0; i < size; i++) {
            if (i != root) MPI_Send(buffer, count, datatype, i, 0, comm);
        }
    } else {
        MPI_Recv(buffer, count, datatype, root, 0, comm, MPI_STATUS_IGNORE);
    }
    return MPI_SUCCESS;
}

int myMPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        for (int i = 0; i < size; i++) {
            if (i == root) {
                memcpy(recvbuf, (char*)sendbuf + i * sendcount * sizeof(int), sendcount * sizeof(int));
            } else {
                MPI_Send((char*)sendbuf + i * sendcount * sizeof(int), sendcount, sendtype, i, 0, comm);
            }
        }
    } else {
        MPI_Recv(recvbuf, recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
    }
    return MPI_SUCCESS;
}

int myMPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        memcpy((char*)recvbuf + rank * recvcount * sizeof(int), sendbuf, sendcount * sizeof(int));
        for (int i = 0; i < size; i++) {
            if (i != root) {
                MPI_Recv((char*)recvbuf + i * recvcount * sizeof(int), recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
            }
        }
    } else {
        MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
    }
    return MPI_SUCCESS;
}

int myMPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    for (int i = 0; i < size; i++) {
        if (i == rank) {
            memcpy((char*)recvbuf + i * recvcount * sizeof(int), (char*)sendbuf + i * sendcount * sizeof(int), sendcount * sizeof(int));
        } else {
            MPI_Send((char*)sendbuf + i * sendcount * sizeof(int), sendcount, sendtype, i, 0, comm);
        }
    }
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Recv((char*)recvbuf + i * recvcount * sizeof(int), recvcount, recvtype, i, 0, comm, MPI_STATUS_IGNORE);
        }
    }
    return MPI_SUCCESS;
}

int myMPI_Reduce(int *sendbuf, int *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        for (int i = 0; i < count; i++) recvbuf[i] = sendbuf[i];

        int *tmp = malloc(count * sizeof(int));
        for (int i = 0; i < size; i++) {
            if (i != root) {
                MPI_Recv(tmp, count, datatype, i, 0, comm, MPI_STATUS_IGNORE);
                for (int j = 0; j < count; j++) {
                    recvbuf[j] += tmp[j]; // Assume MPI_SUM
                }
            }
        }
        free(tmp);
    } else {
        MPI_Send(sendbuf, count, datatype, root, 0, comm);
    }
    return MPI_SUCCESS;
}

int main(int argc, char *argv[]) {
    int rank, size;
    const int root = 0;
    const int n = 1000; // tamanho buffer para testes

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *sendbuf = malloc(n * sizeof(int));
    int *recvbuf = malloc(n * sizeof(int));
    int *alltoall_sendbuf = malloc(n * size * sizeof(int));
    int *alltoall_recvbuf = malloc(n * size * sizeof(int));
    if (!sendbuf || !recvbuf || !alltoall_sendbuf || !alltoall_recvbuf) {
        if (rank == root) printf("Erro alocando memória\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < n; i++) {
        sendbuf[i] = rank + i;
        recvbuf[i] = 0;
    }
    for (int i = 0; i < n * size; i++) {
        alltoall_sendbuf[i] = rank + i;
        alltoall_recvbuf[i] = 0;
    }

    double start, end;
    double mpi_bcast_time, mympi_bcast_time;
    double mpi_scatter_time, mympi_scatter_time;
    double mpi_gather_time, mympi_gather_time;
    double mpi_alltoall_time, mympi_alltoall_time;
    double mpi_reduce_time, mympi_reduce_time;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    MPI_Bcast(sendbuf, n, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mpi_bcast_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    myMPI_Bcast(sendbuf, n, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mympi_bcast_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    MPI_Scatter(sendbuf, n / size, MPI_INT, recvbuf, n / size, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mpi_scatter_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    myMPI_Scatter(sendbuf, n / size, MPI_INT, recvbuf, n / size, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mympi_scatter_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    MPI_Gather(recvbuf, n / size, MPI_INT, sendbuf, n / size, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mpi_gather_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    myMPI_Gather(recvbuf, n / size, MPI_INT, sendbuf, n / size, MPI_INT, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mympi_gather_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    MPI_Alltoall(alltoall_sendbuf, n, MPI_INT, alltoall_recvbuf, n, MPI_INT, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mpi_alltoall_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    myMPI_Alltoall(alltoall_sendbuf, n, MPI_INT, alltoall_recvbuf, n, MPI_INT, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mympi_alltoall_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    MPI_Reduce(sendbuf, recvbuf, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mpi_reduce_time = end - start;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    myMPI_Reduce(sendbuf, recvbuf, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    end = MPI_Wtime();
    mympi_reduce_time = end - start;

    if (rank == root) {
        printf("%-15s | %-15s\n", "Rotina MPI", "Tempo (segundos)");
        printf("-----------------|-----------------\n");
        printf("%-15s | %15.6f\n", "MPI_Bcast", mpi_bcast_time);
        printf("%-15s | %15.6f\n", "myMPI_Bcast", mympi_bcast_time);
        printf("%-15s | %15.6f\n", "MPI_Scatter", mpi_scatter_time);
        printf("%-15s | %15.6f\n", "myMPI_Scatter", mympi_scatter_time);
        printf("%-15s | %15.6f\n", "MPI_Gather", mpi_gather_time);
        printf("%-15s | %15.6f\n", "myMPI_Gather", mympi_gather_time);
        printf("%-15s | %15.6f\n", "MPI_Alltoall", mpi_alltoall_time);
        printf("%-15s | %15.6f\n", "myMPI_Alltoall", mympi_alltoall_time);
        printf("%-15s | %15.6f\n", "MPI_Reduce", mpi_reduce_time);
        printf("%-15s | %15.6f\n", "myMPI_Reduce", mympi_reduce_time);
    }

    free(sendbuf);
    free(recvbuf);
    free(alltoall_sendbuf);
    free(alltoall_recvbuf);

    MPI_Finalize();
    return 0;
}
