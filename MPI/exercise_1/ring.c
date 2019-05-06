#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define FIRST_RANK 0
#define INIT_VALUE 0

int main(int argc, char *argv[]){

    int my_rank;
    int world_size;
    int magic_num;
    int prev_rank;
    int next_rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    prev_rank = (my_rank - 1 > -1? my_rank -1 : world_size - 1);
    next_rank = (my_rank + 1) % world_size;

    if(my_rank == FIRST_RANK){

        magic_num = INIT_VALUE;
	MPI_Send(&magic_num, 1, MPI_INT, next_rank, 1, MPI_COMM_WORLD);
	MPI_Recv(&magic_num, 1, MPI_INT, prev_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("Rank %d received %d from rank %d\n", my_rank, magic_num, prev_rank);

    }else{
	
        MPI_Recv(&magic_num, 1, MPI_INT, prev_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %d received %d from rank %d\n", my_rank, magic_num, prev_rank);
	magic_num = magic_num + 1;
	MPI_Send(&magic_num, 1, MPI_INT, next_rank, 1, MPI_COMM_WORLD);
    }


    MPI_Finalize();
    return 0;

}
