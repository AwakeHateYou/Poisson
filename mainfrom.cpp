#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include <string.h>


#define PI 3.141592
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
	//fout.precision(3);
	for (int j = 0; j < N; j++){
		for (int i = 0; i < N; i++)
			fout << dx*i << ' ' << dx*j << ' ' << mass[i + j*N] << endl;
		fout << endl;
	}
	fout.close();
}

double fillcomplex(double element, double m, double n, double dx){
	return ((((-1.0)/(4.0*N*N)) * element)/(sin(PI*m*dx/N)*sin(PI*m*dx/N) + sin(PI*n*dx/N)*sin(PI*n*dx/N)));
	//return 2*element;
}
double sinmagic(double i, double j){
	return (sin(PI*i/N)*sin(PI*i/N) + sin(PI*j/N)*sin(PI*j/N));
}

int main(){
	double CN = N*((N / 2) + 1);
	double dx = 1.0/N;
	fftw_complex *out, *ncomplex;
	double *U;
	U = new double[N*N];
	initU(U, dx);
	printGNU(U, dx, (char*)"first");
	fftw_plan p, b;
	ncomplex = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	p = fftw_plan_dft_r2c_2d(N, N, U, out, FFTW_ESTIMATE);
	fftw_execute(p);

	double *out_double, *ncomplex_double;
	//for(int k = 0; k < N; k++){
	for(int i = 0; i < N ; i++)
		for(int j = 0; j < N/2+1; j++) {
			out_double = (double*)(&(out[i + j*N]));
			ncomplex_double = (double*)(&(ncomplex[i + j*N]));
			//out_double[0] = fillcomplex(out_double[0], j, i, dx);
			//out_double[1] = fillcomplex(out_double[1], j, i, dx);
			//cout << sinmagic(i, j)<<endl;
			//cout << out_double[0] <<" "<< out_double[1]<< endl;
		    	ncomplex_double[0] = (((-1.0/(4.0*N*N))*out_double[0])/sinmagic(i*dx, j*dx));
			ncomplex_double[1] = (((-1.0/(4.0*N*N))*out_double[1])/sinmagic(i*dx, j*dx));
		}

	b = fftw_plan_dft_c2r_2d(1, N, ncomplex, U, FFTW_ESTIMATE);
	fftw_execute(b);
	fftw_destroy_plan(p);
	fftw_destroy_plan(b);
	fftw_free(out);
	fftw_free(ncomplex);
	printGNU(U, dx, (char*)"second");
	delete[]U;
	return 0;

}
