//
// Poisson 
//
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <cufft.h>


#define PI 3.141592
#define delta 0.1
using namespace std;

#define BLOCK_SIZE 	16			// submatrix size
#define	N			128		// matrix size is N*N

double S(double  x, double y){
	return ((x - 1.0 / 2.0)*(x - 1.0 / 2.0) + (y - 1.0 / 2.0)*(y - 1.0 / 2.0));
}
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
double sinmagic(double i, double j){
	return (sin(PI*i / N)*sin(PI*i / N) + sin(PI*j / N)*sin(PI*j / N));
}
// KERNEL //
__global__ void matMult(cufftDoubleComplex *complex, int n, int dx){
	int   bx = blockIdx.x;		// block index
	int   by = blockIdx.y;
	int   tx = threadIdx.x;		// thread index
	int   ty = threadIdx.y;		

	int   idx = BLOCK_SIZE * bx + tx;
	int	  idy = BLOCK_SIZE * by + ty;

	double *out_double;
	if (idx < (n / 2 + 1) && idy < n){
		out_double = (double*)(&(complex[idx + idy*(n / 2 + 1)]));
		if (idx == 0 && idy == 0)
		{
			out_double[0] = 1;
			out_double[1] = 1;
		}
		else {
			out_double[0] = (((-1.0 / (4.0*n*n))*out_double[0]) /
				(sin(PI*idx / n)*sin(PI*idx / n) + sin(PI*idy / n)*sin(PI*idy / n))) / (n*n);
			out_double[1] = (((-1.0 / (4.0*n*n))*out_double[1]) /
				(sin(PI*idx / n)*sin(PI*idx / n) + sin(PI*idy / n)*sin(PI*idy / n))) / (n*n);
		}
	}
}


// HOST CODE //
int main(int argc, char *  argv[]){
	double dx = 1.0 / N;
	int	numBytesD = N * N * sizeof (cufftDoubleReal);
	int numBytesC = N * (N / 2 + 1) * sizeof (cufftDoubleComplex);

	// allocate host memory
	double *U = new double[N*N];
	initU(U, dx);
	printGNU(U, dx, (char*)"first");
	// allocate device memory
	cufftDoubleReal *Ug;
	cufftDoubleComplex  *complex;
	cufftHandle ahead, backward;

	cudaMalloc((void**)&complex, numBytesC);
	cudaMalloc((void**)&Ug, numBytesD);

	cufftPlan2d(&ahead, N, N, CUFFT_D2Z);
	cufftPlan2d(&backward, N, N, CUFFT_Z2D);



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

	cufftExecD2Z(ahead, Ug, complex);
	matMult <<<blocks, threads >>> (complex, N, dx);
	cufftExecZ2D(backward, complex, Ug);
	

	cudaMemcpy(U, Ug, numBytesD, cudaMemcpyDeviceToHost);
	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gpuTime, start, stop);

	// print the cpu and gpu times
	printf("time spent executing by the GPU: %.2f millseconds\n", gpuTime);
	for (int i = 0; i < N*N; i++)
		U[i] = U[i] / (N*N);

	double shift = U[0];
	for (int i = 0; i < N*N; i++)
		U[i] = U[i] - shift;
	printGNU(U, dx, (char*)"second");
	// release resources
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(Ug);
	cudaFree(complex);
	cufftDestroy(ahead);
	cufftDestroy(backward);

	delete[]U;


	return 0;
}
