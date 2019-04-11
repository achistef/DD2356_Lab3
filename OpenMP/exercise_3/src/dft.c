#include "dft.h"

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int DFT(int idft, double* xr, double* xi, double* Xr_o, double* Xi_o, int N)
{
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for (int k = 0 ; k < N; k++) {
			for (int n = 0; n < N; n++) {
				#pragma omp atomic
				Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
				#pragma omp atomic
				Xi_o[k] += -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
			}
		}
	 	
		// normalize if you are doing IDFT
        	if (idft == -1){
			#pragma omp for
        	        for (int n = 0; n < N; n++) {
                	        Xr_o[n] /= N;
        	                Xi_o[n] /= N;
        	        }
       		}	
	}
	return 1;
}
