#include <stdio.h>
#include <stdlib.h>




void ising(int* G, double* w, int k, int n) {

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

int main() {


	return 0;
}