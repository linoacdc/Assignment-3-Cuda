
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>





__global__ void cudaising(int* G, double* w, int* newG) {

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	double newSpin = 0.0;
	for (int ii = -2; ii <= 2; ii++) {
		for (int jj = -2; jj <= 2; jj++) {

			newSpin += w[(jj + 2) + (ii + 2) * 5] * G[((jj + threadIdx.x + blockDim.x) % blockDim.x) + ((blockIdx.x + ii + blockDim.x) % blockDim.x) * blockDim.x];
		}
	}


	if (newSpin > 0.000001) {
		newG[index] = 1;
	}
	//if newSpin < 0 then the updated spin = -1
	else if (newSpin < -0.000001) {
		newG[index] = -1;
	}
	//if newSpin = 0 then the updated spin = old spin
	else {
		newG[index] = G[index];
	}


}

void cising(int* G, double* w, int k, int n) {

	int* d_G, * d_nG, * temp;
	double* d_w;
	int size = n * n * sizeof(int);
	cudaMalloc((void**)&d_G, size);
	cudaMalloc((void**)&d_w, 5 * 5 * sizeof(double));
	cudaMalloc((void**)&d_nG, size);

	cudaMemcpy(d_G, G, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_w, w, 5 * 5 * sizeof(double), cudaMemcpyHostToDevice);
	//For every timestep
	for (int timestep = 0; timestep < k; timestep++) {

		cudaising <<< n, n >> > (d_G, d_w, d_nG);

		//swap newG with G
		temp = d_nG;
		d_nG = d_G;
		d_G = temp;
		//Now new data are stored in oldG and we can overwrite newG
	}
	cudaMemcpy(G, d_G, size, cudaMemcpyDeviceToHost);
	cudaFree(d_G); cudaFree(d_nG); cudaFree(d_w);
}
