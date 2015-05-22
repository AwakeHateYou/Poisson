#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>

#define N 10
#define delta 0.1
using namespace std;

double S(int x, int y){
	return (pow((x - 1/2), 2) + pow((y - 1/2), 2));
}
void initU(double *mass, double h){
	for (int j = 0; j < N; j++)
		for (int i = 0; i < N; i++)
			if (i == N/2 && j == N/2)
			U[i + N*j] = 0.0;
}
void Print(double *mass){
	ofstream fout;
	fout.open("output.txt");
	fout.precision(3);
	for (int i = 0; i < N*N; i++){
		if ((i % N == 0) && (i != 0))
			fout << std::endl;
		fout << mass[i] << ' ';
	}
	fout.close();
}
void PrintGNU(double *mass){
	ofstream fout;
	fout.open("outputGNU.txt");
	fout.precision(3);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++)
			fout << i << ' ' << j << ' ' << mass[i + j*N] << endl;
		fout << endl;
	}
	fout.close();
}
int main(){
	double h = 1/N;\
	double CN = N*((N/2) +1)
	fftw_complex *out;
	double *U;
	U = new double[N*N];
	initU(U, h);
	fftw_plan p;

	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*CN);
	
	p = fftw_plan_dft_r2c_2d(N*N, CN, U, out, FFTW_ESTIMATE);

	fftw_execute(p);

	fftw_plan b = fftw_plan_dft_c2r_2d(CN, N*N, out, U, FFTW_ESTIMATE);

	fftw_execute(b);

	fftw_destroy_plan(p);
	fftw_destroy_plan(b);
	fftw_free(out);
	PrintGNU(U);
	delete []U;
	return 0;


//	double *U;
//	U = new double[N*N];
//
//	for (int j = 0; j < N; j++){
//		for (int i = 0; i < N; i++)
//			U[i + N*j] = exp(-S(i, j) / (2 * pow(delta, 2)));
//	}
//	//Print(U);
//	PrintGNU(U);
//	delete []U;
//	return 0;
}
