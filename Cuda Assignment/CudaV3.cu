
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>




__global__ void cudaising(int* G, double* w, int* newG, int n, int workperthread) {

	int startingId = threadIdx.x * workperthread;

	//shared w and G in block
	__shared__ double tempW[5 * 5];
	__shared__ int tempG[(517 + 4) * 5];

	//copy necessary elements from G into tempG
	for (int i = -2; i <= 2; i++) {
		for (int j = -2; j <= n + 2; j++) {
			tempG[(j + 2) + (i + 2) * (n + 4)] = G[((j + n) % n) + ((blockIdx.x + i + n) % n) * n];
		}
	}


	//copy using threads
	/*if (threadIdx.x >=25&&threadIdx.x <30) {
		for (int j = -2; j <= n + 2; j++) {
			tempG[(j + 2) + (threadIdx.x-2-25 + 2) * (n + 4)] = G[((j + n) % n) + ((blockIdx.x + threadIdx.x-2-25 + n) % n) * n];

		}
	}
	*/


	//Copy w in tempW


	if (threadIdx.x < 25) {
		tempW[threadIdx.x] = w[threadIdx.x];
	}
	__syncthreads();




	//for every element computed by this thread
	for (int element = 0; element < workperthread; element++) {

		double newSpin = 0.0;

		//for every point in matrix w
		for (int ii = 0; ii < 5; ii++) {
			for (int jj = 0; jj < 5; jj++) {

				//compute new Spin of element
				newSpin += tempW[(jj)+(ii) * 5] * tempG[startingId + element + jj + ii * (n + 4)];

			}
		}
		//global index of element whose spin was just calculated
		int index = startingId + element + blockIdx.x * blockDim.x * workperthread;
		//if newSpin > 0 then the updated spin = 1
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
	__syncthreads();
}



void cising(int* G, double* w, int k, int n) {

	//device variables
	int* d_G, * d_nG, * temp;
	double* d_w;
	int size = n * n * sizeof(int);
	//Allocating space in device for the matrices
	cudaMalloc((void**)&d_G, size);
	cudaMalloc((void**)&d_w, 5 * 5 * sizeof(double));
	cudaMalloc((void**)&d_nG, size);
	//copying the matrices in device
	cudaMemcpy(d_G, G, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_w, w, 5 * 5 * sizeof(double), cudaMemcpyHostToDevice);
	//For every timestep
	for (int timestep = 0; timestep < k; timestep++) {
		//run cudaising with 517 blocks, NOT threads each, NOE elements for each thread
		cudaising << < n, 47 >> > (d_G, d_w, d_nG, n, 11);
		//swap newG with G
		temp = d_nG;
		d_nG = d_G;
		d_G = temp;
		//Now new data are stored in oldG and we can overwrite newG

	}
	//copying final matrix to host
	cudaMemcpy(G, d_G, size, cudaMemcpyDeviceToHost);
	cudaFree(d_G); cudaFree(d_nG); cudaFree(d_w);


}
