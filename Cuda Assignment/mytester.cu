
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


void ising(int* G, double *w, int k, int n) {

	int* oldG = (int*)malloc(n * n * sizeof(int));
	for (int i = 0; i < n * n; i++) {
		oldG[i] = G[i];
	}
	int* newG = (int*)malloc(n * n * sizeof(int));
	double newSpin = 0;
	int* temp;
	//For every timestep
	for (int timestep = 0; timestep < k; timestep++) {

		//for every element
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {

				//get spin of element from the others around it
				for (int ii = -2; ii <= 2; ii++) {
					for (int jj = -2; jj <= 2; jj++) {

						newSpin += w[(jj + 2) + (ii + 2) * 5] * oldG[((jj + j + n) % n) + ((i + ii + n) % n) * n];
					}
				}
				//In newSpin now there is the sum of all spins*weights of all nearby elements

				//if newSpin > 0 then the updated spin = 1
				if (newSpin > 0.000001) {
					newG[j + i * n] = 1;
				}
				//if newSpin < 0 then the updated spin = -1
				else if (newSpin < -0.000001) {
					newG[j + i * n] = -1;
				}
				//if newSpin = 0 then the updated spin = old spin
				else {
					newG[j + i * n] = oldG[j + i * n];
				}
				//reset newSpin
				newSpin = 0;
			}
		}
		//swap newG with G
		temp = newG;
		newG = oldG;
		oldG = temp;
		//Now new data are stored in oldG and we can overwrite newG

	}

	for (int i = 0; i < n * n; i++) {
		G[i] = oldG[i];
	}

}



int main()
{	
	clock_t start, end;
	double cpu_time;
	double w[25] = { 0.004,0.016,0.026 ,0.016 ,0.004 ,0.016,0.071,0.117,0.071,0.016,0.026,0.117,0,0.117, 0.026,0.016, 0.071 , 0.117 ,0.071 ,0.016 ,0.004 ,0.016 ,0.026 ,0.016 ,0.004 };
	int sum = 0;
	int n = 517;
	int k =1;
	int* G = (int*)malloc(n * n * sizeof(int));
	int* cG = (int*)malloc(n * n * sizeof(int));
	
	for (int i = 0; i < n * n; i++) {
		if (rand() % 2 == 0) {
			G[i] = 1;
		}
		else
			G[i] = -1;
	}

	
	memcpy(cG, G, n * n * sizeof(int));

	
	ising(G, w, k, n);
	
	

	//start = clock();
	//cising(cG, w, k, n);
	//end = clock();
	//cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

	for (int i = 0; i < n * n; i++) {
		if (G[i] != cG[i])
			sum++;
	}
	printf("\nErrors: %d\n", sum);
	
	
	//printf("Time with %d moments and %d timesteps: %f seconds\n", n, k, cpu_time);
    return 0;
}