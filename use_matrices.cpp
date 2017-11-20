#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include <ctime>
using namespace std;

int main_linear_solvers()
{
	for (int i=122; i<=200; i+=2){
		// output is like |size of problem = condition number | iters Jacobi | time Jacobi | ...
		// iters GaussSeidel| time GS | iters SOR(1.5) | time SOR(1.5) | ...
		// iters SOR(0.5) | time (SOR0.5)| time LU solve | time QR solve

		////VERSION ITERATIVE
		// cout<<i<<",";
		// Vector x(i);
		// x= randv(x);
		// Matrix A(i,i);
		// A = randm(A);
		// Vector d(i);
		// for (double j=1; j<=i;  j++){
			// d.set_val(int(j),i);
		// }
		// Matrix Lambda(i,i);
		// Lambda= diag(d);
		// A = A +Lambda;

		////VERSION SAME KAPPA k0

		cout<<i<<",";
		Vector x(i);
		x= randv(x);
		Matrix Q(i,i);
		Q = randm(Q);
		Q= qr_q(Q);
		Matrix QT=Q.T();
		Vector d(i);

		double k0 = 2.0;
		for (double j=1; j<=i;  j++){
			double x = k0/j;
			d.set_val(int(j),x);
		}
		Matrix Lambda(i,i);
		Lambda= diag(d);
		Matrix A( i,i);
		A= QT * Lambda;
		A= A*Q;

		////VERSION Diff  KAPPA same size
		// cout<<i<<",";
		// int mmm=50;
		// Vector x(mmm);
		// x= randv(x);
		// Matrix Q(mmm,mmm);
		// Q = randm(Q);
		// Q= qr_q(Q);
		// Matrix QT=Q.T();
		// Vector d(mmm);

		// double step = (i-1.0)/ (mmm-1);
		// for (double j=1; j<=mmm;  j++){
			// double x =1 +step *(j-1);
			// d.set_val(int(j),x);
		// }
		// Matrix Lambda(mmm,mmm);
		// Lambda= diag(d);
		// Matrix A(mmm,mmm);
		// A= QT * Lambda;
		// A= A*Q;
		/////////////////////////////
		Vector b = A*x;
		int count = 0;
		/////////////////////////////////////////////////////
		//clock_t begin = clock();
		//Vector x1 = Jacobi(A,b,count);
		//clock_t end = clock();
		//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//cout<<count<<",";
		//cout<< elapsed_secs<< ",";


		//begin = clock();
		//x1 = SOR(A,b,1.0 , count);
		//end = clock();
		//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//cout<< count<< ","<< elapsed_secs << ",";

		//begin = clock();
		//count = 0;
		//x1 = SOR(A,b,1.5 , count);
	    //end = clock();
	    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//cout<< count<< ","<< elapsed_secs << ",";

		//begin = clock();
		//count = 0;
		//x1 = SOR(A,b,0.5 , count);
		//end = clock();
		//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//cout<< count<< ","<< elapsed_secs << ",";

		clock_t begin = clock();
		count = 0;
		Vector x1 = SD(A,b, count);
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout<< count<<","<< elapsed_secs << ",";

		begin = clock();
		count = 0;
		x1 = CG(A,b, count);
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout<< count<<","<< elapsed_secs << ",";


		begin = clock();
		x1 = LUsolve(A,b);
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		double delx = (x1 - x).norm();
		cout<< elapsed_secs << ","<< delx <<",";

	  begin = clock();
		x1 = QRsolve2(A,b);
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		delx = (x1 - x).norm();
		cout<< elapsed_secs << ","<< delx <<",";

		begin = clock();
		x1 = QRsolve(A,b);
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		delx = (x1 - x).norm();
		cout<< elapsed_secs << ","<< delx;

		cout<<endl;
	}
    return 1;
}

int main()
{
	int i = 20;
	/*
	Matrix G(N,N);
	G = er_graph(G,0.2);
  //cout<<G<<",";
	Matrix D(N,N);
	D = floyd_war(G);
	cout<<D<<",";
	*/
	cout<< endl;
	Vector x(i);
	x= randv(x);
	Matrix Q(i,i);
	Q = randm(Q);
	Q= qr_q(Q);
	Matrix QT=Q.T();
	Vector d(i);

	double k0 = 2.0;
	for (double j=1; j<=i;  j++){
		double x = k0/j;
		d.set_val(int(j),x);
	}
	Matrix Lambda(i,i);
	Lambda= diag(d);
	Matrix A(i,i);
	A= QT * Lambda;
	A= A*Q;

	Vector b = A*x;
	clock_t begin = clock();
	int count = 0;
	Vector x1 = SD(A,b, count);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< count<<","<< elapsed_secs << ",";

	begin = clock();
	count = 0;
	x1 = momentum(A,b, count,0.5,0.2);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< count<<","<< elapsed_secs << ",";


	cout<< endl;
	return 1;
}
