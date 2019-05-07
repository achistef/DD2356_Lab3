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
    int iter;
    for (iter = 0; iter < flip / num_ranks; iter++)
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

    int total_count;
    MPI_Reduce(local_count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0){
		// Estimate Pi and display the result
		*answer = ((double)total_count / (double)flip) * 4.0;
	}
}
