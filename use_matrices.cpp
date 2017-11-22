#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include <ctime>
using namespace std;
#include <cmath>
#include <fstream>

int main()
{
	ofstream outfile ("same_kappa_new.txt");
	int i = 1;
	for (int j=5; j<=60; j+=1){
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
    i = pow(1.08,j)+10;
		outfile<<i<<",";
  //
	  double sd_co =0;
		double sd_ti =0;
		double cg_co =0;
		double cg_ti =0;
		double cgp_co =0;
		double cgp_ti =0;
		double lu_ti =0;
		double lu_dx =0;
		double qr_ti =0;
		double qr_dx =0;
		double qre_ti =0;
		double qre_dx =0;

	  for (int k=0;k<=10;k++){
			Vector x(i);
			x = randv(x);
			Matrix Q(i,i);
			Q = randm(Q);
			Q = qr_q(Q);
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
			A = QT * Lambda;
			A = A*Q;

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
	    Matrix D = diag(A);
	    //// code starts here /////////////////////////////////////
			clock_t begin = clock();
			count = 0;
			Vector x1 = SD(A,b, count);
			clock_t end = clock();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      sd_ti += elapsed_secs;
			sd_co += count;

			begin = clock();
			count = 0;
			x1 = CG(A,b, count);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			cg_ti += elapsed_secs;
			cg_co += count;

			begin = clock();
			count = 0;
			x1 = CG_pre(A,b, count);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			cgp_ti += elapsed_secs;
			cgp_co += count;

			//////////////////// non-iterative solvers
			begin = clock();
			x1 = LUsolve(A,b);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			double delx = (x1 - x).norm();
			lu_ti += elapsed_secs;
			lu_dx += delx;

		  begin = clock();
			x1 = QRsolve2(A,b);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			delx = (x1 - x).norm();
			qr_ti += elapsed_secs;
			qr_dx += delx;

			begin = clock();
			x1 = QRsolve(A,b);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			delx = (x1 - x).norm();
			qre_ti += elapsed_secs;
			qre_dx += delx;
		}

		outfile<< sd_co/10<<","<< sd_ti/10 << ",";
		outfile<< cg_co/10<<","<< cg_ti/10 << ",";
		outfile<< cgp_co/10<<","<< cgp_ti/10 << ",";
		outfile<< lu_ti/10<<","<< lu_dx/10 << ",";
		outfile<< qr_ti/10<<","<< qr_dx/10 << ",";
		outfile<< qre_ti/10<<","<< qre_dx/10;
		outfile<<endl;
		outfile.flush();
	}
	outfile.close();
  return 1;
}

int main_graph()
{
	int N = 600;
	Matrix G(N,N);
	G = er_graph(G,0.7);
	G = G - diag(G);
  //cout<<G<<",";
	Matrix D(N,N);
	//D = floyd_war(G);
	//cout<<D<<",";
	Vector e(N);
	e = evec_c(G);
	cout << e<<endl;
	// Vector x(i);
	// x = randv(x);
	// Matrix Q(i,i);
	// Q = randm(Q);
	// Q = qr_q(Q);
	// Matrix QT=Q.T();
	// Vector d(i);
  //
	// double k0 = 2.0;
	// for (double j=1; j<=i;  j++){
	// 	double x = k0/j;
	// 	d.set_val(int(j),x);
	// }
	// Matrix Lambda(i,i);
	// Lambda= diag(d);
	// Matrix A(i,i);
	// A= QT * Lambda;
	// A= A*Q;
  //
	// Vector b = A*x;
	// clock_t begin = clock();
	// int count = 0;
	// Vector x1 = SD(A,b, count);
	// clock_t end = clock();
	// double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	// cout<< count<<","<< elapsed_secs << ",";
  //
	// begin = clock();
	// count = 0;
	// x1 = momentum(A,b, count,0.5,0.2);
	// end = clock();
	// elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	// cout<< count<<","<< elapsed_secs << ",";

  cout<<
	cout<< endl;
	return 1;
}
