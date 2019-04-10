#include "sum.h"

void omp_sum(double *sum_ret)
{
	double sum = 0;
	#pragma omp parallel shared(sum)
	{
		int id = omp_get_thread_num();
		int thread_count = omp_get_num_threads();
		for(int i = id;  i < size; i+=thread_count){
			sum+= x[i];
		}
	}
	*sum_ret = sum;
}

void omp_critical_sum(double *sum_ret)
{
	double sum = 0;
        #pragma omp parallel shared(sum)
        {
                int id = omp_get_thread_num();
                int thread_count = omp_get_num_threads();
                for(int i = id;  i < size; i+=thread_count){
                	#pragma omp critical
			sum+= x[i];
                }
        }
        *sum_ret = sum;
}

void omp_atomic_sum(double *sum_ret)
{
 	double sum = 0;
        #pragma omp parallel shared(sum)
        {
                int id = omp_get_thread_num();
                int thread_count = omp_get_num_threads();
                for(int i = id;  i < size; i+=thread_count){
                	#pragma omp atomic
                        sum+= x[i];
                }
        }
        *sum_ret = sum;
}

void omp_local_sum(double *sum_ret)
{
 	// initialize an array that will store partial sums
        int thread_count = omp_get_max_threads();
        double *sum = (double*)calloc(thread_count, sizeof(double));
        #pragma omp parallel
        {
                // each thread calculates only a portion of sums
                int id = omp_get_thread_num();
                for(int i = id;  i < size; i+=thread_count){
                        sum[id]+= x[i];
                }
        }

        // serially sum up the array
        double acc = 0;
        for(int i = 0; i < thread_count; i++){
                acc += sum[i];
        }
        *sum_ret = acc;
}

void omp_padded_sum(double *sum_ret)
{
 	// initialize an array that will store partial sums
        int thread_count = omp_get_max_threads();
	// add 56 Bytes of padding
        double *sum = (double*)calloc(thread_count, 8*sizeof(double));
        #pragma omp parallel
        {
                // each thread calculates only a portion of sums
                int id = omp_get_thread_num();
                for(int i = id;  i < size; i+=thread_count){
                        sum[8*id]+= x[i];
                }
        }

        // serially sum up the array
        double acc = 0;
        for(int i = 0; i < 8 * thread_count; i++){
                acc += sum[i];
        }
        *sum_ret = acc;
}

void omp_private_sum(double *sum_ret)
{
	*sum_ret = 0;
	double sum;
	// initialize an array that will store partial sums
        int thread_count = omp_get_max_threads();
        #pragma omp parallel private(sum)
        {
		sum = 0;
                // each thread calculates only a portion of sums
                int id = omp_get_thread_num();
                for(int i = id;  i < size; i+=thread_count){
                        sum+= x[i];
                }

		#pragma omp critical
		*sum_ret += sum;
        }
}

void omp_reduction_sum(double *sum_ret)
{
        double sum = 0;
        // initialize an array that will store partial sums
        #pragma omp parallel for reduction(+:sum)
        for(int i = 0;  i < size; i++) sum+= x[i];
	*sum_ret = sum;

}
