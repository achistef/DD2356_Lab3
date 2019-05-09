#include "block_matmul.h"

struct Config {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile;

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *B_tmp, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2];
	int B_dims[2];
	int C_dims[2];
	int matrix_size;

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;


void init_matmul(char *A_file, char *B_file, char *outfile)
{

	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);

	/* Copy output file name to configuration */
	config.outfile = (char*) malloc(sizeof(char) * (strlen(outfile)+1));
   	strcpy(config.outfile, outfile);


	/* Get matrix size header */
	if(config.world_rank == 0){
        	MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
        	MPI_File_read(config.A_file, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        	MPI_File_close(&config.A_file);
		
		MPI_File_open(MPI_COMM_SELF, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);
                MPI_File_read(config.B_file, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_close(&config.B_file);
	}


	/* Broadcast global matrix sizes */
	MPI_Bcast(config.A_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(config.B_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
	config.dim[0] = config.dim[1] = sqrt(config.world_size);


	/* Verify dim of A and B matches for matul and both are square*/
	if(config.A_dims[0] != config.A_dims[1]){
                MPI_Finalize();
        }

	if(config.B_dims[0] != config.B_dims[1]){
                MPI_Finalize();
        }

	if(config.A_dims[1] != config.B_dims[0]){
		MPI_Finalize();
	}
	
	if(config.A_dims[0] % config.dim[0] != 0 || config.A_dims[1] % config.dim[1] != 0 || config.B_dims[0] % config.dim[0] != 0 || config.B_dims[1] % config.dim[1] != 0){
		MPI_Finalize();
	}


	/* Create Cart communicator for NxN processes */
	int period[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period , 1, &config.grid_comm);

	MPI_Cart_coords(config.grid_comm, config.world_rank, 2, config.coords);
        //printf("Rank %d coordinates are %d %d\n", config.world_rank, config.coords[0], config.coords[1]);

	MPI_Cart_rank(config.grid_comm, config.coords, &config.grid_rank);
        //printf("The processor at position (%d, %d) has rank %d\n", config.coords[0], config.coords[1], config.grid_rank);
	

	/* Sub div cart communicator to N row communicator */
	int selected_dim[2] = {0, 1};
	MPI_Cart_sub(config.grid_comm, selected_dim, &config.row_comm);
	MPI_Comm_rank(config.row_comm, &config.row_rank);
	MPI_Comm_size(config.row_comm, &config.row_size);
	printf("Rank %d: Row comm rank %d, world size %d\n", config.world_rank, config.row_rank, config.row_size);


	/* Sub div cart communicator to N col communicator */
	selected_dim[0] = 1;
	selected_dim[1] = 0;
        MPI_Cart_sub(config.grid_comm, selected_dim, &config.col_comm);
	MPI_Comm_rank(config.col_comm, &config.col_rank);
        MPI_Comm_size(config.col_comm, &config.col_size);
	printf("Rank %d: Col comm rank %d, world size %d\n", config.world_rank, config.col_rank, config.col_size);
	

	/* Setup sizes of full matrices */
	config.C_dims[0] = config.A_dims[0];
	config.C_dims[1] = config.B_dims[1];
	config.matrix_size = config.C_dims[0] * config.C_dims[1];

	/* Setup sizes of local matrix tiles */
	config.local_dims[0] = config.local_dims[1] = config.C_dims[0] / config.dim[0];
	config.local_size = config.local_dims[0] * config.local_dims[1];


	/* Create subarray datatype for local matrix tile */
	int start_from[2] = {config.row_rank * config.local_dims[0], config.col_rank * config.local_dims[1]};
	//printf("Rank %d is responsible for %d %d\n",config.world_rank, start_from[0], start_from[1]);
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, start_from, MPI_ORDER_C, MPI_DOUBLE, &config.block);
	MPI_Type_commit(&config.block);


	/* Create data array to load actual block matrix data */
	config.A     = (double*) malloc(config.local_size* sizeof(double));
	config.A_tmp = (double*) malloc(config.local_size* sizeof(double));
	config.B     = (double*) malloc(config.local_size* sizeof(double));
	config.B_tmp = (double*) malloc(config.local_size* sizeof(double));
	config.C     = (double*) calloc(config.local_size, sizeof(double));

	/* Set fileview of process to respective matrix block */
	MPI_File_open(MPI_COMM_WORLD, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
	MPI_File_open(MPI_COMM_WORLD, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);

	MPI_File_set_view(config.A_file, 2*sizeof(int), MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, 2*sizeof(int), MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

	/* Collective read blocks from files */
	MPI_File_read_all(config.A_file, config.A, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, config.B, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);

	// print first four numbers of each process
	//printf("Rank %d read A : %f %f %f %f...\n", config.world_rank, config.A[0], config.A[1], config.A[2], config.A[3]);
	//printf("Rank %d read B : %f %f %f %f...\n", config.world_rank, config.B[0], config.B[1], config.B[2], config.B[3]);	
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */
	if(config.world_rank == 0){
		MPI_File_open(MPI_COMM_SELF, config.outfile, MPI_MODE_WRONLY | MPI_MODE_CREATE , MPI_INFO_NULL, &config.C_file);
		MPI_File_write(config.C_file, config.C_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_close(&config.C_file);
	}


	/* Set fileview of process to respective matrix block with header offset */
	MPI_File_open(MPI_COMM_WORLD, config.outfile, MPI_MODE_WRONLY, MPI_INFO_NULL, &config.C_file);
	MPI_File_set_view(config.C_file, 2*sizeof(int), MPI_DOUBLE, config.block, "native", MPI_INFO_NULL);

	/* Collective write and close file */
	MPI_File_write_all(config.C_file, config.C, config.local_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&config.C_file);

	/* Cleanup */
}

void multiply_matrices()
{
 	int i, j, k ;
	int dim = config.local_dims[0];
  	for (i = 0 ; i < dim ; i++) {
    		for (k = 0 ; k < dim ; k++) {
      			for (j = 0 ; j < dim ; j++) {
				config.C[(i * dim) + j] += config.A_tmp[(i * dim) + k] * config.B[(k * dim) + j];
      			}
    		}
  	}
}

void array_copy(double *a, double*b, int count){
	int i;
	for(i = 0; i < count; i++){
		a[i] = b[i];
	}
}

void compute_fox()
{
	//printf("Rank %d B: %f %f %f %f...\n", config.world_rank, config.B[0], config.B[1], config.B[2], config.B[3]);

	/* Compute source and target for verticle shift of B blocks */
	int above = config.col_rank -1;
	if(above <0) above +=config.dim[0];
	int below = (config.col_rank +1) % config.dim[0];

	int r; //round 
	for (r = 0; r < config.dim[0]; r++) {
		
		// calculate which col rank will bcast
		int bcast_col = (config.col_rank + r) % config.dim[0];

		if(bcast_col == config.row_rank){
			array_copy(config.A_tmp, config.A, config.local_size);
			//if (config.world_rank == 0) 
			//	printf("Rank %d A tmp before bcast: %f %f %f %f...\n", config.world_rank, config.A_tmp[0], config.A_tmp[1], config.A_tmp[2], config.A_tmp[3]);
			//printf("Round %d Rank %d will bcast!\n", r, config.world_rank);
		}
		
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		MPI_Bcast(config.A_tmp,  config.local_size,  MPI_DOUBLE,  bcast_col,  config.row_comm);
		//if (config.world_rank == 0)
			//printf("Rank %d A tmp after bcast : %f %f %f %f...\n", config.world_rank, config.A_tmp[0], config.A_tmp[1], config.A_tmp[2], config.A_tmp[3]);

		/* dgemm with blocks */
		multiply_matrices();
		//if (config.world_rank == 0)
		//printf("Rank %d output: %f %f %f %f...\n", config.world_rank, config.C[0], config.C[1], config.C[2], config.C[3]);	
		//if (config.world_rank == 0)
                //        printf("Rank %d C mult : %f %f %f %f...\n", config.world_rank, config.C[0], config.C[1], config.C[2], config.C[3]);	
		
		/* Shfting block B upwards and receive from process below */
		array_copy(config.B_tmp, config.B, config.local_size);
		MPI_Request request;
    		MPI_Status status;
		MPI_Isend(config.B_tmp, config.local_size, MPI_DOUBLE, above, r+1, config.col_comm, &request);
		MPI_Recv(config.B, config.local_size, MPI_DOUBLE, below, r+1, config.col_comm, MPI_STATUS_IGNORE);
		MPI_Wait(&request, &status);
		////printf("Rank %d output: %f %f %f %f...\n", config.world_rank, config.C[0], config.C[1], config.C[2], config.C[3]);
		//printf("Rank %d B: %f %f %f %f...\n", config.world_rank, config.B[0], config.B[1], config.B[2], config.B[3]);
		//MPI_Barrier(MPI_COMM_WORLD);
		///printf("End\n");
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rank %d C: %f %f %f %f...\n", config.world_rank, config.C[0], config.C[1], config.C[2], config.C[3]);

}
