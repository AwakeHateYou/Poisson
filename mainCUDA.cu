//
// Perform "naive" square matrix multiplication
//
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include <string.h>


#define PI 3.141592
#define delta 0.1
using namespace std;

#define BLOCK_SIZE 	16			// submatrix size
#define	N			128		// matrix size is N*N

void initU(double *mass, double dx){
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++)
		mass[i + N*j] = ((S(dx*i, dx*j) - 2.0 * delta*delta) / (delta*delta*delta*delta))
		*exp(-S(dx*i, dx*j) / 2.0*delta*delta);
}
void printGNU(double *mass, double dx, char* filename){
	ofstream fout;
	fout.open(filename);
	//fout.precision(3);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++)
			fout << dx*i << ' ' << dx*j << ' ' << mass[i + j*N] << endl;
		fout << endl;
	}
	fout.close();
}

// KERNEL //
__global__ void matMult(cufftDoubleReal *U, int n){
	int   bx = blockIdx.x;		// block index
	int   by = blockIdx.y;
	int   tx = threadIdx.x;		// thread index
	int   ty = threadIdx.y;
	float sum = 0.0f;			// computed subelement
	int	  ia = n * BLOCK_SIZE * by + n * ty;	// a [i][0]
	int   ib = BLOCK_SIZE * bx + tx;

	// Multiply the two matrices together;
	for (int k = 0; k < n; k++)
		sum += a[ia + k] * b[ib + k*n];

	// Write the block sub-matrix to global memory;
	// each thread writes one element
	int ic = n * BLOCK_SIZE * by + BLOCK_SIZE * bx;

	c[ic + n * ty + tx] = sum;
}


// HOST CODE //
int main(int argc, char *  argv[]){
	double dx = 1.0 / N;
	int	numBytesD = N * N * sizeof (cufftDoubleReal);
	int numBytesC = N * (N / 2 + 1) * sizeof (cufftDoubleComplex);

	// allocate host memory
	double *U = new double[N*N];
	initU(U, dx);

	// allocate device memory
	cufftDoubleReal *Ug;
	cufftDoubleComplex  *out, *ncomplex;
	cudaMalloc((void**)&out, numBytesC);
	cudaMalloc((void**)&ncomplex, numBytesC);
	cudaMalloc((void**)&Ug, numBytesD);


	// set kernel launch configuration
	dim3 threads(BLOCK_SIZE, BLOCK_SIZE);
	dim3 blocks(N / threads.x, N / threads.y);

	// create cuda event handles
	cudaEvent_t start, stop;
	float gpuTime = 0.0f;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// asynchronously issue work to the GPU (all to stream 0)
	cudaEventRecord(start, 0);
	cudaMemcpy(Ug, U, numBytesD, cudaMemcpyHostToDevice);

	matMult << <blocks, threads >> > (Ug, N);

	cudaMemcpy(c, cdev, numBytes, cudaMemcpyDeviceToHost);
	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gpuTime, start, stop);

	// print the cpu and gpu times
	printf("time spent executing by the GPU: %.2f millseconds\n", gpuTime);

	// release resources
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(out);
	cudaFree(Ug);
	cudaFree(ncomplex);

	delete []U;


	return 0;
}
