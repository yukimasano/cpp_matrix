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
	ofstream outfile;
	int i = 1;

	for (int j=15; j<=50; j+=1){

		/* format of output is
		# format    0 N , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time,
		#           5 CG_pre count 6 CG_pre time
		#           7 LU time, 8 LU delx, 9 QR time, 10 QR delx
		#						11 full QR time 12 full QR delx
		*/

		double sd_co =0.0;
		double sd_ti =0.0;
		double cg_co =0.0;
		double cg_ti =0.0;
		double cgp_co =0.0;
		double cgp_ti =0.0;
		double lu_ti =0.0;
		double lu_dx =0.0;
		double qr_ti =0.0;
		double qr_dx =0.0;
		double qre_ti =0.0;
		double qre_dx =0.0;
    i = pow(1.09,j)+8;
		outfile<<i<<",";



		// for loop to average measurements
	  for (int k=0;k<=10 ;k++){
			////// VERSION SAME KAPPA DIFFERENT SIZE N
			clock_t begin = clock();
			Vector x(i);
			x = randv(x);
			Matrix Q(i,i);
			Q = rand_basis_gs(Q);
			Matrix QT=Q.T();
			Matrix A(i,i);
			A = Q;
			double k0 = 2.0;
			for (int jj=1; jj<=i;  jj++){
				double x = ((k0-1)/float(i-1))*float(jj) + 1.;
				for (int h =1; h<=i; h++){
					A.set_val(jj,h, A(jj,h)*x);
				}
			}
			Vector b = A*x;
			clock_t end = clock();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			cout<<"dim = "<<i<<" time for prep: "<<elapsed_secs<<endl;

			int count = 0;
	    //// code starts here /////////////////////////////////////
			begin = clock();
			count = 0;
			Vector x1 = SD(A,b, count);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
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

			////////////////// non-iterative solvers
			begin = clock();
			x1 = LUsolve(A,b);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			double delx = (x1 - x).norm();
			lu_ti += elapsed_secs;
			lu_dx += delx;

			if (i<200){
				//cout<<"QR2"<<flush<<endl;
			  begin = clock();
				x1 = QRsolve2(A,b);
				end = clock();
				elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				delx = (x1 - x).norm();
				qr_ti += elapsed_secs;
				qr_dx += delx;
			}
			else{
				qr_ti = 0;
				qr_dx = 0;
			}
			//cout<<"QRs"<<flush<<endl;
      if (i<100){
				begin = clock();
				x1 = QRsolve(A,b);
				end = clock();
				elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				delx = (x1 - x).norm();
				qre_ti += elapsed_secs;
				qre_dx += delx;
			}
			else{
				qre_ti = 0;
				qre_dx = 0;
			}
		}
		//cout<<"done"<<flush<<endl;
		outfile.open("varyN_kappa2.txt", fstream::app);
		outfile<< sd_co/10<<","<< sd_ti/10 << ",";
		outfile<< cg_co/10<<","<< cg_ti/10 << ",";
		outfile<< cgp_co/10<<","<< cgp_ti/10 << ",";
		outfile<< lu_ti/10<<","<< lu_dx/10 << ",";
		outfile<< qr_ti/10<<","<< qr_dx/10 << ",";
		outfile<< qre_ti/10<<","<< qre_dx/10<<",";
		outfile<<endl;
		outfile.flush();
		outfile.close();
	}

  return 1;
}
