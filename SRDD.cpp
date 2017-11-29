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

	for (int j=16; j<=110; j+=1){
		/* format of output is
		# format    0 kappa=size , 1 SD Count, 2 SD Time, 3 CG Count, 4 CG time,
		#           5 CG_pre count 6 CG_pre time
		#           7 LU time, 8 LU delx, 9 QR time, 10 QR delx
		#						11 full QR time 12 full QR delx
		#						13 Jacobi count 14 Jacobi Time
		#						15 SOR1 "-" 16
		#						17 SOR1.5 "-" 18
		#						19 SOR0.5 "-" 20
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
		double jac_co =0.0;
		double jac_ti =0.0;
		double sor1_co =0.0;
		double sor1_ti =0.0;
		double sor15_co =0.0;
		double sor15_ti =0.0;
		double sor05_co =0.0;
		double sor05_ti =0.0;
    i = pow(1.09,j)+8;




		// for loop to average measurements
	  for (int k=0;k<=10 ;k++){
			//VERSION ITERATIVE
			clock_t begin = clock();
			Vector x(i);
			x= randv(x);
			Matrix A(i,i);
			A = randm(A);
			for (int o=1; o<=i;  o++){
				A.set_val(o,o,o+i);
			}
			Vector b = A*x;
			Matrix AA = A;

			clock_t end = clock();
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			cout<<"dim = "<<i<<" time for prep: "<<elapsed_secs<<endl;
			// /////////////////////////////
			int count = 0;
	    //// code starts here /////////////////////////////////////
			Vector x1(i);
			begin = clock();
			count = 0;
			x1 = SD(A,b, count);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      sd_ti += elapsed_secs;
			sd_co += count;
			cout<<"SD"<<flush<<endl;

			// else{
			// 	sd_ti = 0;
			// 	sd_co = 0;
			// }

			// begin = clock();
			// count = 0;
			// x1 = CG(AA,b2, count);
			// end = clock();
			// elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			// cg_ti += elapsed_secs;
			// cg_co += count;

			// cout<<"CG"<<flush<<endl;
      //
			// begin = clock();
			// count = 0;
			// x1 = CG_pre(AA,b2, count);
			// end = clock();
			// elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			// cgp_ti += elapsed_secs;
			// cgp_co += count;


			////////////////// non-iterative solvers

			begin = clock();
			x1 = LUsolve(A,b);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			double delx = (x1 - x).norm();
			lu_ti += elapsed_secs;
			lu_dx += delx;
      //
			// if (i<150){
			// 	//cout<<"QR2"<<flush<<endl;
			// 	cout<<"QR2"<<flush<<endl;
      //
			//   begin = clock();
			// 	x1 = QRsolve2(A,b);
			// 	end = clock();
			// 	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			// 	delx = (x1 - x).norm();
			// 	qr_ti += elapsed_secs;
			// 	qr_dx += delx;
			// }
			// else{
			// 	qr_ti = 0;
			// 	qr_dx = 0;
			// }
      //
			// cout<<"QR"<<flush<<endl;
      // if (i<200){
			// 	begin = clock();
			// 	x1 = QRsolve(A,b);
			// 	end = clock();
			// 	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			// 	delx = (x1 - x).norm();
			// 	qre_ti += elapsed_secs;
			// 	qre_dx += delx;
			// }
			// else{
			// 	qre_ti = 0;
			// 	qre_dx = 0;
			// }
			cout<<"Jac"<<flush<<endl;

			count = 0;
			begin = clock();
			x1 = Jacobi(A,b,count);
			end = clock();
			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			jac_ti += elapsed_secs;
			jac_co += count;

			cout<<"SOR1"<<flush<<endl;

			begin = clock();
			if (i<140){
				x1 = SOR(A,b,1.0 , count);
				end = clock();
				elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				sor1_ti += elapsed_secs;
				sor1_co += count;
				cout<<"SOR15"<<flush<<endl;

				begin = clock();
				count = 0;
				x1 = SOR(A,b,1.5 , count);
		    end = clock();
		    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				sor15_ti += elapsed_secs;
				sor15_co += count;
			}
			if (i<100){
				cout<<"SOR5"<<flush<<endl;
				begin = clock();
				count = 0;
				x1 = SOR(A,b,0.5 , count);
				end = clock();
				elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				sor05_ti += elapsed_secs;
				sor05_co += count;
		  }

		}
		//cout<<"done"<<flush<<endl;
		outfile.open("SRDD.txt", fstream::app);
		outfile<<i<<",";
		outfile<< sd_co/10<<","<< sd_ti/10 << ",";
		outfile<< cg_co/10<<","<< cg_ti/10 << ",";
		outfile<< cgp_co/10<<","<< cgp_ti/10 << ",";
		outfile<< lu_ti/10<<","<< lu_dx/10 << ",";
		outfile<< qr_ti/10<<","<< qr_dx/10 << ",";
		outfile<< qre_ti/10<<","<< qre_dx/10<<",";
		outfile<< jac_co/10<<","<< jac_ti/10 << ",";
		outfile<< sor1_co/10<<","<< sor1_ti/10 << ",";
		outfile<< sor15_co/10<<","<< sor15_ti/10 << ",";
		outfile<< sor05_co/10<<","<< sor05_ti/10;
		outfile<<endl;
		outfile.flush();
		outfile.close();
	}

  return 1;
}
