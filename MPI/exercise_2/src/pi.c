#include "pi.h"
#include <math.h>
#include <time.h>


void init_pi(int set_seed, char *outfile)
{
	if (filename != NULL) {
		free(filename);
		filename = NULL;
	}

	if (outfile != NULL) {
		filename = (char*)calloc(sizeof(char), strlen(outfile)+1);
		memcpy(filename, outfile, strlen(outfile));
		filename[strlen(outfile)] = 0;
	}
	seed = set_seed;
}

void cleanup_pi()
{
	if (filename != NULL)
		free(filename);
}

void compute_pi(int flip, int *local_count, double *answer)
{
    double x, y, z;
    int rank;
    int num_ranks;
    *local_count = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    srand(seed * rank); // Important: Multiply SEED by "rank" when you introduce MPI!

    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < flip / num_ranks; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            (*local_count) ++;
        }
    }

    if (rank == 0){
    	long local_counts[num_ranks];
    	MPI_Request requests[num_ranks];

		for (int src = 1; src < num_ranks; src++) 
			MPI_Irecv(&local_counts[src - 1], 1, MPI_INT, src, 0, MPI_COMM_WORLD, &requests[src - 1]);
		
		MPI_Waitall(num_ranks - 1, requests, MPI_STATUS_IGNORE);

		int total_count = *local_count;
		for (int i = 1; i < num_ranks; i ++) 
			total_count += local_counts[i - 1];

		// Estimate Pi and display the result
		*answer = ((double)total_count / (double)flip) * 4.0;
	}
	else {
		MPI_Request request;
		MPI_Isend(local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
	}
}
