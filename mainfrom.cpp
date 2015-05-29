#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include <string.h>


#define PI 3.14159265358979323846
#define N 100
#define delta 0.1
using namespace std;


double S(double  x, double y){
	return ((x - 1.0 / 2.0)*(x - 1.0 / 2.0) + (y - 1.0 / 2.0)*(y - 1.0 / 2.0));
}
void initU(double *mass, double dx){
	for (int j = 0; j < N; j++)
	for (int i = 0; i < N; i++)
		mass[i + N*j] = ((S(dx*i, dx*j) - 2.0 * delta*delta )/ (delta*delta*delta*delta))
		*exp(-S(dx*i, dx*j) / 2.0*delta*delta);
		//mass[i + N*j] = exp(-S(dx*i, dx*j) / (2 * pow(delta, 2)));
//	for (int j = 0; j < N; j++)
//	for (int i = 0; i < N; i++){
//		if (i == 0)
//			mass[i + j*N] = mass[N - 1 + j*N];
//		if (j == 0)
//			mass[i + j*N] = mass[i + (N - 1)*N];
//	}


}
void print(double *mass, double dx){
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
void printGNU(double *mass, double dx, char* filename){
	ofstream fout;
	fout.open(filename);
	fout.precision(3);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++)
			fout << dx*i << ' ' << dx*j << ' ' << mass[i + j*N] << endl;
		fout << endl;
	}
	fout.close();
}

//fftw_complex VVV(double m, double n, fftw_complex* R) {
//	return (((-1/4*N*N)*R(m,n))/(math.sin(PI*m/N)*math.sin(PI*m/N) + math.sin(PI*n/N)*math.sin(PI*n/N)));
//}

int main(){
	double CN = N*((N / 2) + 1);
	double dx = 1.0/N;
	fftw_complex *out, *Uc;
	double *U;
	U = new double[N*N];
	initU(U, dx);
	printGNU(U, dx, (char*)"before.txt");
	fftw_plan p, b;
	Uc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	p = fftw_plan_dft_r2c_2d(N, N, U, out, FFTW_ESTIMATE);
	fftw_execute(p);

	double *D, *B;

	for(int i = 0; i < N ; i++)
		for(int j = 0; j < N/2 + 1; j++) {
			D = (double*)(&(out[i + j*N]));
			B = (double*)(&(Uc[i + j*N]));
			//cout << D[0] <<" "<< D[1]<< endl;
		    B[0] = (((-1.0/(4.0*N*N))*D[0])/(sin(PI*i*dx/N)*sin(PI*i*dx/N) + sin(PI*j*dx/N)*sin(PI*j*dx/N)));
			B[1] = (((-1.0/(4.0*N*N))*D[1])/(sin(PI*i*dx/N)*sin(PI*i*dx/N) + sin(PI*j*dx/N)*sin(PI*j*dx/N)));
		}

//	for(int i = 0; i < N/2 + 1; i++)
//		for(int j = 0; j < N; j++) {
//			B = (double*)(&(Uc[i + j*N]));
//			cout << B[0] << " " << B[1] << endl;
//		}

	b = fftw_plan_dft_c2r_2d(N, N, Uc, U, FFTW_ESTIMATE);
	fftw_execute(b);
	fftw_destroy_plan(p);
	fftw_destroy_plan(b);
	fftw_free(out);
	fftw_free(Uc);
	printGNU(U, dx, (char*)"after.txt");
	delete[]U;
	return 0;

}
