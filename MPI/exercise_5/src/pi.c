#include "pi.h"
#include <math.h>
#include <time.h>

void write_in_file();
int size_of_array();
void add_new_value();

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


    int prev_rank = (rank - 1 > -1? rank -1 : num_ranks - 1);
    int next_rank = (rank + 1) % num_ranks;
    int token = 0;

    if(rank == 0){

		MPI_Send(&token, 1, MPI_INT, next_rank, 1, MPI_COMM_WORLD);
		MPI_Recv(&token, 1, MPI_INT, prev_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		write_in_file(rank, (double)(*local_count) / ((double)flip / (double)num_ranks));

    }else{
	
        MPI_Recv(&token, 1, MPI_INT, prev_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        write_in_file(rank, (double)(*local_count) / ((double)flip / (double)num_ranks));
		token = token + 1;
		MPI_Send(&token, 1, MPI_INT, next_rank, 1, MPI_COMM_WORLD);
    }

    if (rank == 0){   	
    	// rank 0 is the last process that accesses the file. 
    	// We are sure that all other processes have finished writing.
		MPI_File fh;
		MPI_Status status;
		
		MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

		char ch;
		char buffer[100];
		buffer[0] = '\0';
		double total_count = 0;

		while(1) {

			MPI_File_read(fh, &ch, 1, MPI_CHAR, &status);

			int count = 0;
			MPI_Get_count(&status, MPI_CHAR, &count);
			if(count < 1) 
				break;		//EOF
			else 
				if(ch != '\n')
					strcat(buffer, &ch);
				else {
					add_new_value(buffer, &total_count, num_ranks);
					buffer[0] = '\0';
				}
		}
	
		*answer = (double)total_count * 4.0;

		char pi[100];
		sprintf (pi, "pi = %f\n", *answer);
		MPI_File_write(fh, pi, size_of_array(pi), MPI_CHAR, MPI_STATUS_IGNORE);
		
		MPI_File_close(&fh);
	}
}

void add_new_value(char *buffer, double *total_count, int num_ranks) {
	char *value_token;
	double value;
	char sep[2] = " ";
	strtok(buffer, sep);
	value_token = strtok(NULL, sep);
	value = atof(value_token);
	*total_count += value / (double)num_ranks;
}

void write_in_file(int rank, double value){
	char result[100]; 
	sprintf (result, "%d %f\n", rank, value);
	MPI_File fh;
	MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
	MPI_File_write(fh, result, size_of_array(result), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
}

int size_of_array(char* pointer){
	int i = 0;
	while(pointer[i] != '\0') i++;
	return i;
}
