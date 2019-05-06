#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]){


    MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    printf("hello there\n");


    MPI_Finalize();
    return 0;

}