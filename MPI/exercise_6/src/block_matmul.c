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
	MPI_Comm grid_comm;				//set
	MPI_Comm row_comm;				//set
	MPI_Comm col_comm;				//set

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];				//set
	int world_rank, world_size, grid_rank;		//set
	int row_rank, row_size, col_rank, col_size;	//set

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
	int selected_dim[2] = {1, 0};
	MPI_Cart_sub(config.grid_comm, selected_dim, &config.row_comm);
	MPI_Comm_rank(config.row_comm, &config.row_rank);
	MPI_Comm_size(config.row_comm, &config.row_size);
	//printf("Rank %d: Row comm rank %d, world size %d\n", config.world_rank, config.row_rank, config.row_size);


	/* Sub div cart communicator to N col communicator */
	selected_dim[0] = 0;
	selected_dim[1] = 1;
        MPI_Cart_sub(config.grid_comm, selected_dim, &config.col_comm);
	MPI_Comm_rank(config.col_comm, &config.col_rank);
        MPI_Comm_size(config.col_comm, &config.col_size);
	//printf("Rank %d: Col comm rank %d, world size %d\n", config.world_rank, config.col_rank, config.col_size);
	

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
	config.A = malloc(config.local_size*sizeof(double));
	config.A_tmp = malloc(config.local_size*sizeof(double));
	config.B = malloc(config.local_size*sizeof(double));
	config.C = malloc(config.local_size*sizeof(double));

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

	//PRINT
	printf("Rank %d read: %f %f %f %f...\n", config.world_rank, config.A[0], config.A[1], config.A[2], config.A[3]);
	
}

void cleanup_matmul()
{
	printf(" ");
	/* Rank zero writes header specifying dim of result matrix C */

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void compute_fox()
{
	printf(" ");

	/* Compute source and target for verticle shift of B blocks */
	int i;
	for (i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */

		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */

	}
}
