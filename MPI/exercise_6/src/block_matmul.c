#include "block_matmul.h"

struct Config {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile; 					//set

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;		//set 1,2
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2]; 					//set
	int B_dims[2]; 					//set
	int C_dims[2];
	int matrix_size;

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;

void init_matmul(char *A_file, char *B_file, char *outfile)
{
	printf("init matmul...\n");
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
	printf("rank %d: my dim is %d", config.world_rank, config.dim[0]);


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


	/* Create Cart communicator for NxN processes */
	int period[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, period , 1, &config.grid_comm);
	

	/* Sub div cart communicator to N row communicator */

	/* Sub div cart communicator to N col communicator */

	/* Setup sizes of full matrices */

	/* Setup sizes of local matrix tiles */

	/* Create subarray datatype for local matrix tile */

	/* Create data array to load actual block matrix data */

	/* Set fileview of process to respective matrix block */

	/* Collective read blocks from files */

	/* Close data source files */
}

void cleanup_matmul()
{
	printf("cleanup...\n");
	/* Rank zero writes header specifying dim of result matrix C */

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void compute_fox()
{
	printf("compute fox...\n");

	/* Compute source and target for verticle shift of B blocks */
	int i;
	for (i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */

		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */

	}
}
