
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>



__global__ void cudaising(int* G, double* w, int* newG, int workperthread, int n) {

	//global id of first element whose spin is to be computed by this thread
	int startingId = (threadIdx.x + blockIdx.x * blockDim.x) * workperthread;
	int element = 0;
	//for each element which will be calculated by this thread
	for (int element = 0; element < workperthread; element++) {
		//calculate spin of element with global id = starting id + element
		double newSpin = 0.0;
		for (int ii = -2; ii <= 2; ii++) {
			for (int jj = -2; jj <= 2; jj++) {

				newSpin += w[(jj + 2) + (ii + 2) * 5] * G[((jj + threadIdx.x * workperthread + n + element) % n) + ((blockIdx.x + ii + n) % n) * n];
			}
		}
		int index = startingId + element;
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
}



void cising(int* G, double* w, int k, int n) {
	//
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
		//cudaising with 517 blocks and 47 threads each and 11 moments calculated per thread
		cudaising << < n, 47 >> > (d_G, d_w, d_nG, 11, n);

		//swap newG with G
		temp = d_nG;
		d_nG = d_G;
		d_G = temp;
		//Now new data are stored in oldG and we can overwrite newG
	}
	//copy updated matrix to host
	cudaMemcpy(G, d_G, size, cudaMemcpyDeviceToHost);
	cudaFree(d_G); cudaFree(d_nG); cudaFree(d_w);


}
