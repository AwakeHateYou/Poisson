#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>



#define N 100
#define delta 0.1
#define dx 0.01
using namespace std;


double S(double  x, double y){
	return (pow((x - 1 / 2.0), 2) + pow((y - 1 / 2.0), 2));
}
void initU(double *mass){
	for (int j = 1; j < N; j++)
	for (int i = 1; i < N; i++)
		mass[i + N*j] = ((S(dx*i, dx*j) - 2 * pow(delta, 2)) / pow(delta, 4))
		*exp(-S(dx*i, dx*j) / (2 * pow(delta, 2)));
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++){
		if (i == 0)
			mass[i + j*N] = mass[N - 1 + j*N];
		if (j == 0)
			mass[i + j*N] = mass[i + (N - 1)*N];
	}


}
void print(double *mass){
	ofstream fout;
	fout.open("output.txt");
	//fout.precision(3);
	for (int i = 0; i < N*N; i++){
		if ((i % N == 0) && (i != 0))
			fout << std::endl;
		fout << mass[i] << ' ';
	}
	fout.close();
}
void printGNU(double *mass){
	ofstream fout;
	fout.open("outputGNU.txt");
	fout.precision(3);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++)
			fout << dx*i << ' ' << dx*j << ' ' << mass[i + j*N] << endl;
		fout << endl;
	}
	fout.close();
}
int main(){
	double CN = N*((N / 2) + 1);
	fftw_complex *out;
	double *U;
	U = new double[N*N];
	initU(U);
	printGNU(U);
	fftw_plan p, b;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	p = fftw_plan_dft_r2c_2d(N - 1, N - 1, U, out, FFTW_ESTIMATE);
	fftw_execute(p);
	b = fftw_plan_dft_c2r_2d(N - 1, N - 1, out, U, FFTW_ESTIMATE);
	fftw_execute(b);
	fftw_destroy_plan(p);
	fftw_destroy_plan(b);
	fftw_free(out);
	printGNU(U);
	delete[]U;
	return 0;

}