#include <iostream>
#include <cmath>
#include <fstream>
#define N 10
#define delta 0.1
using namespace std;

double S(int x, int y){
	return pow((x - 1/2), 2) + pow((y - 1/2), 2);
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
	fout.open("output.txt");
	fout.precision(3);
	for (int i = 0; i < N*N; i++){
		if ((i % N == 0) && (i != 0))
		fout << mass[i] << ' ';
	}
	fout.close();
}
int main(){
	double *U;
	U = new double[N*N];
	for (int j = 0; j < N; j++)
		for (int i = 0; i < N; i++)
			U[i + N*j] = exp(-S(i, j) / 2 * delta);
		Print(U);

	return 0;
}