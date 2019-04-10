#include "sum.h"

void omp_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_critical_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_atomic_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_local_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_padded_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_private_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}

void omp_reduction_sum(double *sum_ret)
{
	double sum = 0;

	for (int i = 0; i < size; i++){
		sum += x[i];
	}

	*sum_ret = sum;
}
