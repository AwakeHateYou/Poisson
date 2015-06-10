#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>
#include <string.h>


#define PI 3.141592
#define N 16
#define delta 0.1
using namespace std;


double S(double  x, double y){
	return ((x - 1.0 / 2.0)*(x - 1.0 / 2.0) + (y - 1.0 / 2.0)*(y - 1.0 / 2.0));
}
void initU(double *mass, double dx){
	for (int j = 0; j < N; j++)
		for (int i = 0; i < N; i++)
			mass[i + N*j] = ((S(dx*i, dx*j) - 2.0 * pow(delta, 2.0) )/ (pow(delta, 4.0)))
			*exp(-S(dx*i, dx*j) / 2.0*pow(delta, 2.0));
	
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
		for (int i = 0; i < N; i++) {
//			cout << "M = " << mass[i+j*N] << endl;
			fout << dx*i << ' ' << dx*j << ' ' << mass[i + j*N] << endl;
		}
		fout << endl;
	}
	fout.close();
}

double sinmagic(double i, double j){        
	return (sin(PI*i/N)*sin(PI*i/N) + sin(PI*j/N)*sin(PI*j/N));
}

int main(){
	int CN = N*((N / 2) + 1);
	double dx = 1.0/N;
	double constant = (-1.0/(4.0*N*N));
	fftw_complex *out;
	double *U;
	U = new double[N*N];
	initU(U, dx);
	printGNU(U, dx, (char*)"first.txt");
	fftw_plan p, b;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*CN);
	p = fftw_plan_dft_r2c_2d(N, N, U, out, FFTW_ESTIMATE);
	fftw_execute(p);

	double *out_double;
	for(int i = 0; i < N ; i++)
		for(int j = 0; j < (N/2+1); j++) {
			out_double = (double*)(&(out[i*(N/2+1) + j]));	
		  	if (i==0 && j==0) {	
				out_double[0] = 1;
				out_double[1] = 1;
		  	}
		  	else{
				out_double[0] = (constant*(out_double[0])/sinmagic(i, j))/(N*N);
				out_double[1] = (constant*(out_double[1])/sinmagic(i, j))/(N*N);
			
		  	}
			}


	b = fftw_plan_dft_c2r_2d(N, N, out, U, FFTW_ESTIMATE);
	fftw_execute(b);
	
	
	for (int i = 0; i < N*N; i++)
	    U[i] = U[i]/(N*N);
	
	double shift  = U[0];
	for (int i = 0; i < N*N; i++)
	    U[i] = U[i]-shift;
	
	
	
	printGNU(U, dx, (char*)"second.txt");
	fftw_destroy_plan(p);
	fftw_destroy_plan(b);
	fftw_free(out);
	delete[]U;
	return 0;

}
