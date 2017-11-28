#include <iostream>
#include "Matrix.hpp"
#include <tuple>
#include <cassert>
#include <math.h>
#include <stdlib.h>
// constructor that creates matrix of given size with
// double precision entries all initially set to zero
Matrix::Matrix(int rows, int cols)
{
  mData = new double* [rows];
  mShape[0] = rows;
  mShape[1] = cols;
	for (int i=0; i<rows;  i++){
		mData[i] = new double[cols];
		for (int j=0; j<cols; j++){
			mData[i][j] = 0.0;
  	}
  }
}


// destructor - deletes pointer
Matrix::~Matrix()
{
  for (int i=0; i<mShape[0];  i++){
		delete[] mData[i];
  }
  delete mData;
}

double& Matrix::operator()(int i, int j)
{

  if (i < 1 || j< 1)
    {
      throw Exception("out of range",
		  "accessing vector through () - index too small");
    }
  else if (i > mShape[0] || j >mShape[1] )
    {
      throw Exception("length mismatch",
		  "accessing vector through () - index too high");
    }
  return mData[i-1][j-1];
}

Matrix Matrix::T() const
{
	Matrix C(mShape[1],mShape[0]);
	for (int i=0; i<mShape[1];  i++){
		for (int j=0; j<mShape[0]; j++){
			C.mData[i][j] = mData[j][i];
		}
	}
	return C;
}
int nrows(const Matrix& C)
{
	return C.mShape[0];
}

int ncols(const Matrix& C)
{
	return C.mShape[1];
}

Vector Matrix::get_col(int k)
{
	//std::cout<<"before all"<<std::endl<< std::flush;
	Vector v(mShape[0]);
  //	std::cout<<"after vec"<<std::endl<< std::flush;
	for (int i =1 ; i<=mShape[0]; i++){
		//std::cout<<"pre-definition"<<std::endl<< std::flush;
		//std::cout<< "this is the val  "<< mData[i][k]<< "this is the size  "<<mShape[0] << std::endl<< std::flush;
		v.set_val(i, mData[i-1][k-1]);
		//std::cout<<"in for loop"<<std::endl<< std::flush;
	}
	//std::cout<<"before return"<<std::endl<< std::flush;
	return v;
}

Vector Matrix::get_row(int k)
{
	Vector v(mShape[0]);
	for (int i =1 ; i<=mShape[1]; i++){
		v.set_val(i, mData[k-1][i-1]);
	}
	return v;
}


Matrix operator+(const Matrix& A, const Matrix& B)
{
	assert(A.mShape[1] == B.mShape[1] && A.mShape[0] == B.mShape[0]);
	Matrix W(A.mShape[0], A.mShape[1]);
  	W.mData=new double* [W.mShape[0]];
  	for (int i=0; i<W.mShape[0];  i++){
  		W.mData[i] = new double[W.mShape[1]];
  		for (int j=0; j<W.mShape[1]; j++){
  			W.mData[i][j] = A.mData[i][j]+ B.mData[i][j];
  		}
  	}
  return W;
}

Matrix operator+(const Matrix& A, const double k)
{
	Matrix W(A.mShape[0], A.mShape[1]);
  	W.mData=new double* [W.mShape[0]];
  	for (int i=0; i<W.mShape[0];  i++){
  		W.mData[i] = new double[W.mShape[1]];
  		for (int j=0; j<W.mShape[1]; j++){
  			W.mData[i][j] = A.mData[i][j]+ k;
  		}
  	}
  return W;
}

Matrix operator-(const Matrix& A, const Matrix& B)
{
	assert(A.mShape[1] == B.mShape[1] && A.mShape[0] == B.mShape[0]);
	Matrix W(A.mShape[0], A.mShape[1]);
  	W.mData=new double* [W.mShape[0]];
  	for (int i=0; i<W.mShape[0];  i++){
  		W.mData[i] = new double[W.mShape[1]];
  		for (int j=0; j<W.mShape[1]; j++){
  			W.mData[i][j] = A.mData[i][j]- B.mData[i][j];
  		}
  	}
  return W;
}


// definition of vector operator =
Matrix& Matrix::operator=(const Matrix& A)
{
  if (A.mShape[0] != mShape[0] || A.mShape[1] != mShape[1])
    {
      throw Exception("size mismatch",
		  "matrix assignment operator - matrices have different lengths");
    }
  else
    {
	  for (int i=0; i<mShape[0];  i++){
	  		mData[i] = new double[mShape[1]];
	  		for (int j=0; j<mShape[1]; j++){
	  			mData[i][j] = A.mData[i][j];
	  		}
	  }
  }
  return *this;
}


Matrix operator-(const Matrix& C)
{
	Matrix W(C.mShape[0],C.mShape[1]);
	W.mData=new double* [W.mShape[0]];
	for (int i=0; i<W.mShape[0];  i++){
		W.mData[i] = new double[W.mShape[1]];
		for (int j=0; j<W.mShape[1]; j++){
			W.mData[i][j] = -C.mData[i][j];
		}
	}
	return W;
}

Matrix operator*(Matrix& A, Matrix& B){
	assert(A.mShape[1] == B.mShape[0]);
	Matrix C(A.mShape[0],B.mShape[1]);
	//std::cout<<C<<std::endl;
  Matrix AT = A.T();
	for (int i=0; i<A.mShape[0];  i++){
		for (int j=0; j<B.mShape[1]; j++){
		//	C.mData[i][j] = s_prod(A.get_col(i+1), B.T().get_col(j+1), A.mShape[0]);
			C.set_val(i+1,j+1 , (AT.get_col(i+1) * B.get_col(j+1)));
		}
	}
	return C;
}

Vector operator*(Matrix& A, Vector& b){
	assert(nrows(A) == length(b));
	int k = length(b);
	Vector v(k);
	for (int i=1; i<=k;  i++){
		v.set_val(i, A.get_row(i)*b);
	}
	return v;
}

Matrix operator*(const Matrix& A, const double& a)
{
	Matrix C(A.mShape[0],A.mShape[1]);
	for (int i=0; i<A.mShape[0];  i++){
		for (int j=0; j<A.mShape[1]; j++){
			C.mData[i][j] = a*A.mData[i][j];
		}
	}
	return C;
}

Matrix operator/(const double& a,const Matrix& A)
{
	Matrix C(A.mShape[0],A.mShape[1]);
	for (int i=0; i<A.mShape[0];  i++){
		for (int j=0; j<A.mShape[1]; j++){
      if (A.mData[i][j]!=0){
			     C.mData[i][j] = a/A.mData[i][j];
      }
		}
	}
	return C;
}

Matrix operator*(const double& a, const Matrix& A)
{
	Matrix C(A.mShape[0],A.mShape[1]);
	for (int i=0; i<A.mShape[0];  i++){
		for (int j=0; j<A.mShape[1]; j++){
			C.mData[i][j] = a*A.mData[i][j];
		}
	}
	return C;
}

std::ostream& operator<<(std::ostream& output, const Matrix& C) {
  for (int i=0; i<C.mShape[0]; i++)
    {
	  output << "(";
	  for (int j=0; j<C.mShape[1]; j++)
	      {
		  output <<  C.mData[i][j];
		  if (j<C.mShape[1]-1)
			  output  << ",\t";
		  else
			  output  << "";
	      }
	  output <<")\n";
    }
  return output;  // for multiple << operators.
}

Matrix diag(const Matrix& A){
	Matrix C(A.mShape[0],A.mShape[1]);
	for (int i=0; i<A.mShape[0];  i++){
		C.mData[i][i] = A.mData[i][i];
	}
	return C;
}

Vector Matrix::diag(){
  Vector v(mShape[0]);
	if (mShape[0] <= mShape[1]){
		for (int i=0; i<mShape[0];  i++){
			v.set_val( i+1 , mData[i][i]);
		}
	}
	else {
		Vector v(mShape[1]);
		for (int i=0; i<mShape[1];  i++){
				v.set_val(i+1, mData[i][i]);
			}
		}
	return v;
}


void Matrix::set_val(int i, int j, double x){
	mData[i-1][j-1] = x;
}
//
// void Matrix::set_col(int i, Vector& x){
// 	for (int k=0;k<mShape[0];k++){
//     mData[k][i-1] = x(k);
//   }
// }

Matrix diag(Vector& v){
	int k =length(v);
	Matrix C(k,k);
	for (int i=1; i<=k; i++){
		C.set_val(i,i, v(i));
	}
	return C;
}


Matrix::Matrix(const Matrix& C){
	mShape[0]=C.mShape[0];
	mData=new double* [C.mShape[0]];
	mShape[1]=C.mShape[1];
	for (int i=0; i<C.mShape[0];  i++){
		mData[i] = new double[C.mShape[1]];
			for (int j=0; j<C.mShape[1]; j++){
				mData[i][j] =  C.mData[i][j];
			}
	}
}


Matrix randm(Matrix& A){
	int m = nrows(A);
	//srand( (unsigned)time( NULL ) );
	for (int i=1; i<=m ; i++){
		for (int j=1; j<=m ; j++){
			A.set_val(i,j, (float)rand()/RAND_MAX);
		}
	}
	return A;
}

Matrix rand_basis_gs(Matrix& A){
  // using modified Gram-Schmidt
  assert(nrows(A)==ncols(A));
  int i = nrows(A);
  A = randm(A);
  Vector ve(i);
  for (int k=1;k<=i; k++){
    if (k==1){ve=A.get_col(k);}
    else{
      ve = A.get_col(k);
      // do the Gram-Schmidt
      for (int m=1; m<k;m++){
        ve = ve - (A.get_col(m)*A.get_col(k))* A.get_col(m);
      }
    }
    ve = ve/ve.norm();
    //A.set_col(k,ve);
    for (int j=1;j<=i; j++){
     A.set_val(j,k, ve(j));
    }
  }
return A;
}

Matrix rand_basis_angle(Matrix& A){
  //using spherical coordinates in N-d
  assert(nrows(A)==ncols(A));
  int m = nrows(A);
  A = randm(A);
  Vector d = A.diag();
  A.set_val(1,1, A(1,1)*2.0);
  A = A*M_PI;
  Matrix B = A;
  for (int k=1;k<=m;k++){
    Vector vei = A.get_col(k);
    for (int j=1;j<=m;j++){
      if (j<k){
        vei.set_val(j, -cos(B(j,j))); //orthogonal to previous
      }
      else{
  	    vei.set_val(j,cos(vei(j))); // any direction in that orthogonal space
      }
    }
    for (int j=1;j<=m;j++){
  	   A.set_val(j,k, vei(j)/vei.norm());
    }
  }
return A;
}


Matrix er_graph(Matrix& A,float p){
  assert(nrows(A)==ncols(A));
	int s = nrows(A);
	//srand( (unsigned)time( NULL ) );
	for (int i=1; i<=s ; i++){
		for (int j=1; j<=i ; j++){
			A.set_val(i,j, p>(float)rand()/RAND_MAX);
      A.set_val(j,i, p>(float)rand()/RAND_MAX);
		}
	}
	return A;
}

Matrix floyd_war(Matrix& A){
  assert(nrows(A)==ncols(A));
  int s = nrows(A);
  Matrix B(s,s);
  B= B + 999999999;
	for (int i=1; i<=s; i++){
    B.set_val(i,i,0);
  }
  for (int k=1; k<=s; k++){
		for (int i=1; i<=k; i++){
      if (A(k,i)!=0){
        B.set_val(k,i, A(k,i));
        B.set_val(i,k, A(k,i));
      }
    }
  }
  for (int k=1; k<=s; k++){
		for (int i=1; i<=s; i++){
      for (int j=1; j<=s; j++){
        if (B(i,j) > B(i,k) + B(k,j)){
          B.set_val(i,j,B(i,k)+B(k,j));
        }
      }
		}
	}
	return B;
}
/*
Vector betweenness_c(Matrix& D){
  assert(nrows(D)==ncols(D));
  int s = nrows(D);
  Vector v;
  for (int k=1; k<=nrows(D); k++){
    double a;
    for (int j=1; j<=nrows(D); j++){
      if ((j!=k) && (D(k,j)!=999999999)){
        a+=
      }
    }
  }
	return v;
}
*/
// LINEAR SOLVERS



Vector evec_c(Matrix& A){
	assert(nrows(A) == ncols(A));
	// INITIALISATION
	Vector x(nrows(A));
  Vector xold(nrows(A));
  xold  = randv(xold);
  x  = randv(x);
  float eps = 1.0;
  float lam;
  while (eps>1e-15){
    x = A*x;
    lam = (xold*x)/xold.norm();
    x = x/x.norm();
    eps=(x-xold).norm();
    xold=x;
  }
  signed int s = (lam>0) - (lam<0);
  x= x*s;
  return x;
}

Vector LUsolve(Matrix& A, Vector& v){
	assert(nrows(A) == ncols(A));
	// we want a square system

	// INITIALISATION
	Matrix U=A;
	int m =nrows(A);
	Vector temp=ones(m);
	Matrix L = diag(temp);
	//LU FACTORISATION
	for (int k=1; k<=m-1; k++){
		for(int j = k+1; j<=m; j++){
			L.set_val(j,k,   (U(j,k)/U(k,k)));
			for (int l=k; l<=m; l++){
				U.set_val(j,l, (U(j,l) - U(k,l)*L(j,k)));
			}
		}
	}
	//std::cout<<L<<std::endl;
	Vector y(m);
	// FORWARD SOLVE
	for (int k=0; k<m; k++){
		y.set_val(k+1, v(k+1));
		for(int j = 0; j<k; j++){
			y.set_val(k+1, y(k+1)-(L(k+1,j+1))*y(j+1));
		}
		y.set_val(k+1, y(k+1)/L(k+1,k+1));
		}
	//std::cout<<y<<std::endl;
	Vector x(m);
	// BACKSUBSTITUTION
	for (int k=m-1; k>=0; k--){
		x.set_val(k+1, y(k+1));
		for (int j = k+1; j<m; j++){
			x.set_val(k+1, x(k+1)- U(k+1,j+1)*x(j+1));
		}
		x.set_val(k+1, x(k+1)/U(k+1,k+1));
	}
	//std::cout<<U<<std::endl;
	return x;
}
Vector QRsolve(Matrix& A, Vector& v){
	// QR via Householder
	assert(nrows(A) == ncols(A));
	int m = ncols(A);
	Matrix R = A;
	Vector temp=ones(m);
	Matrix Q = diag(temp);
    for (int j = 0; j<m; j++){
    	Matrix x(m-j,1);
    	double len=0.0;
    	for (int l=j; l<m; l++){
    		x.set_val(l+1-j,1 ,R(l+1,j+1));
    		len+=R(l+1,j+1) * R(l+1,j+1);
    	}

    	len=sqrt(len);

    	int sign = (x(1,1) >= 0) - (x(1,1) < 0);
    	Matrix e(m-j,1);
    	e.set_val(1,1,1.0);

    	Matrix u = (sign*len*e + x);

    	// normalise u
    	len=0.0;
    	for (int i = 1; i<=m-j; i++){
    		len+=u(i,1)*u(i,1);
    	}
      //std::cout<<j<<u<<std::endl<<std::flush;
      if (len!=0){
    	   u=u* (1./sqrt(len));
      }
      //   	std::cout<<"here comes  u " << std::endl<< j<<"  "<< u << std::endl<<std::flush;

    	// SUBMULTIPLICATION    START /////////////////////////////
    	Matrix sub(m-j, m-j);
    	Matrix subQ(m,m-j);
    	for (int i = j; i<m; i++){
    	//	std::cout<<i<<std::flush<<std::endl;
    		for (int h = j; h<m; h++){
    			sub.set_val(i+1 - j,h+1 - j, R(i+1,h+1));
    		}
    	}

    	for (int i = 0; i<m; i++){
			for (int h = j; h<m; h++){
				subQ.set_val(i+1,h+1-j,Q(i+1,h+1));
			}
		}

    	Matrix w=(u.T());
    	Matrix z= w*sub;
      //  	std::cout<<"here is sub pre   "<<std::endl<< u*z<<std::flush<< std::endl;
    	sub = sub -2.0*(u*z);

    	Vector temp2= ones(m-j);
    	Matrix subeye= diag(temp2);
    	Matrix temp3= subeye -2.0*(u*w);
    	subQ= subQ*(temp3);

      //    	std::cout<<"here is sub post   "<<std::endl<< sub<<std::flush<< std::endl;
    	for (int i = j; i<m; i++){
    		for (int h = j; h<m; h++){
    			R.set_val(i+1,h+1, sub(i+1- j,h+1- j));
    		}
    	}

    	for (int i = 0; i<m; i++){
			for (int h = j; h<m; h++){
				Q.set_val(i+1, h+1, subQ(i+1,h+1-j));
			}
      //std::cout<<j<<"end of for"<<std::endl<<std::flush;
		}

    	// SUBMULTIPLICATION  END  ////////////////////////////

    }
    //    std::cout<< "here is R "<< std::endl<< R<<std::endl<<std::flush;
    //    std::cout<< "here is Q "<< std::endl<< Q<<std::endl<<std::flush;
    //    std::cout<< "here is QR "<< std::endl<< Q*R<<std::endl<<std::flush;

    // Invert Q multiply b with it.
    Matrix QT=Q.T();
    Vector v2=QT* v;
    //BACKWARD SOLVE FOR X
	Vector x(m);
	for (int k=m-1; k>=0; k--){
		x.set_val(k+1, v2(k+1));
		for (int j = k+1; j<m; j++){
			x.set_val(k+1, x(k+1)- R(k+1,j+1)*x(j+1));
			}
		x.set_val(k+1, x(k+1)/R(k+1,k+1));
		}
    //	std::cout<< A*x -v <<std::endl;
	return x;
}
Vector QRsolve2(Matrix& A, Vector& bb){
	assert(nrows(A) == ncols(A));
	int m = ncols(A);
  Matrix A2 = A;
  Vector b=bb;
	Matrix V(m-1,m-1);
	for (int k=0; k<m ;k++){
		Vector xx(m-k);
    //std::cout<<k<<xx<<std::flush<<std::endl;
		for (int j=1;j<=m-k; j++){
			xx.set_val(j, A(j+k,k+1));
		}
		//std::cout<<xx<<std::endl<<std::flush;
		Vector vk(m-k);

		double len =xx.norm();
		vk.set_val(1,len);
    //std::cout<<k<<vk<<std::endl<<std::flush;
    //std::cout<<k<<xx.norm()<<std::endl<<std::flush;


		vk = vk+ xx;
    if (vk.norm()!=0){
		    vk = vk / vk.norm();
    }
		Matrix A_sub(m-k,m-k);
    	for (int i = k; i<m; i++){
    		for (int h = k; h<m; h++){
    			A_sub.set_val(i-k+1,h-k+1, A(i+1,h+1));
    		}
    	}
    //std::cout<<k<<A<<std::endl<<std::flush;
    	//std::cout<<A_sub<<std::endl<<std::flush;
		Matrix vk2(m-k,1);
		for (int jj=1 ; jj<=m-k;jj++){
			vk2.set_val(jj,1, vk(jj));
		}
		Matrix vk2T= vk2.T();
		Matrix vk3 = vk2T*A_sub;
		//std::cout<<vk2<<std::endl<<std::flush;
		Matrix A_sub2 = vk2*vk3;


		A_sub2 = 2*A_sub2;
		A_sub2 = A_sub - A_sub2;
		//std::cout<<A<<std::endl<<std::flush;

		for (int j = k; j<m; j++){
			for (int h=k; h<m; h++){
				A2.set_val(j+1,h+1,A_sub2(j+1-k, h+1-k) );
			}
		}
		Vector bk2(m-k);
		for (int jj=k ; jj<m;jj++){
			bk2.set_val(jj+1-k,b(jj+1));
		}
		Vector bk3  =  vk;
		double temp = (vk*bk2);
		//std::cout<<temp<<std::endl<<std::flush;
		bk3 = temp*vk;
		bk3 = 2*bk3;
		bk3 = bk2 - bk3;
		for (int jj=k ; jj<m;jj++){
			b.set_val(jj+1, bk3(jj-k+1));
		}
		for (int j=1; j<=m-k-1; j++){
			V.set_val(j,k+1,vk(j));
		}
	}
	//std::cout<<b<<std::endl<<std::flush;
	Vector x(m);
	/// BACKSOLving////////
	for (int k=m-1; k>=0; k--){
		x.set_val(k+1, b(k+1));
		for (int j = k+1; j<m; j++){
			x.set_val(k+1, x(k+1)- A2(k+1,j+1)*x(j+1));
			}
		x.set_val(k+1, x(k+1)/A2(k+1,k+1));
		}
	return x;
}

Vector Jacobi(Matrix& A, Vector& v, int& count){
	assert(nrows(A) == ncols(A));
	int m =nrows(A);
	// we want a square system

	// INITIALISATION
	Matrix T=A;
	for (int k=1; k<=m; k++){
		T.set_val(k,k, (1./A(k,k)));
	}
	Vector temp = ones(m);
	Matrix eye = diag(temp);
  Vector c = T*v;
	T =  (eye - T);
	double tol = 1.e-14*v.norm();
	Vector x=c;
	//std::cout<<" before the while loop "<<(A*x -c).norm() << std::endl<<std::flush;
  Vector r = (A*x -v);
	while (r.norm() > tol){
		count+=1;
		x = T*x + c;
    r = (A*x -v);
	//	std::cout<<(A*x -v).norm() << std::endl<<std::flush;
	}
  //	std::cout<<count<<std::endl<<std::flush;

	return x;
}
Vector SOR(Matrix& A, Vector& v, double w, int& count){
	assert(nrows(A) == ncols(A));
	int m =nrows(A);
	// we want a square system

	// INITIALISATION
	Matrix D=diag(A);
	for (int k=1; k<=m; k++){
		if (D(k,k)!= 0) {
		D.set_val(k,k, (1./D(k,k)));
		}
	}
	Matrix U(m,m);
	Matrix L(m,m);
	for (int i=1; i<=m; i++){
		for (int j=1; j<=m; j++){
			if (i>j) {
				L.set_val(i,j, A(i,j));
			}
			if (i<j) {
				U.set_val(i,j, A(i,j));
			}
		}
	}
	Vector temp = ones(m);
	Matrix Lbar= D*L;
	Matrix LL= (diag(temp) + w*Lbar);
	// T = (1+wLbar)**-1 ((1-w)11 - wU)
	// T = LL ((1-w)11 - wU)
	Matrix Ubar=D*U;
	Matrix M = diag(temp);
	M= M*(1.-w) -  w*Ubar;
	Vector c =  (D*v) *w;
	double tol = 1.e-14*v.norm();
	Vector x=c;
	Vector rh(m);
  Vector r = c;
	while (r.norm() >tol){
		count+=1;
		rh = M*x + c;
		//std::cout<<x<<std::endl<<std::flush;
		//forward solve
		for (int k=0; k<m; k++){
			x.set_val(k+1, rh(k+1));
			for(int j = 0; j<k; j++){
				x.set_val(k+1, x(k+1)-(LL(k+1,j+1))*x(j+1));
			}
			x.set_val(k+1, x(k+1)/LL(k+1,k+1));
		}
    r = (A*x -v);
	}
  //	std::cout<<count<<std::endl<<std::flush;
	return x;
}

Matrix qr_q(Matrix& A){
	int m = ncols(A);
	Matrix R= A;
	Vector temp=ones(m);
	Matrix Q = diag(temp);
    for (int j = 0; j<m; j++){
    	Matrix x(m-j,1);
    	double len=0.0;
    	for (int l=j; l<m; l++){
    		x.set_val(l+1-j,1 ,R(l+1,j+1));
    		len+=R(l+1,j+1) * R(l+1,j+1);
    	}

    	len=sqrt(len);

    	int sign = (x(1,1) >= 0) - (x(1,1) < 0);
    	Matrix e(m-j,1);
    	e.set_val(1,1,1.0);

    	Matrix u = (sign*len*e + x);

    	// normalise u
    	len=0.0;
    	for (int i = 1; i<=m-j; i++){
    		len+=u(i,1)*u(i,1);
    	}
      if (len!=0){u= (1./len)*u;}
      //   	std::cout<<"here comes  u " << std::endl<< j<<"  "<< u << std::endl<<std::flush;

    	// SUBMULTIPLICATION    START /////////////////////////////
    	Matrix sub(m-j, m-j);
    	Matrix subQ(m,m-j);
    	for (int i = j; i<m; i++){
    		for (int h = j; h<m; h++){
    			sub.set_val(i+1 - j,h+1 - j, R(i+1,h+1));
    		}
    	}

    	for (int i = 0; i<m; i++){
			for (int h = j; h<m; h++){
				subQ.set_val(i+1,h+1-j,Q(i+1,h+1));
			}
		}

    	Matrix w=(u.T());
    	Matrix z= w*sub;
    	sub = sub -2.0*(u*z);

    	Vector temp2= ones(m-j);
    	Matrix subeye= diag(temp2);
    	Matrix temp3= subeye -2.0*(u*w);
    	subQ= subQ*(temp3);
    	for (int i = j; i<m; i++){
    		for (int h = j; h<m; h++){
    			R.set_val(i+1,h+1, sub(i+1- j,h+1- j));
    		}
    	}
      std::cout<<R<<std::flush<<std::endl;
    	for (int i = 0; i<m; i++){
  			for (int h = j; h<m; h++){
  				Q.set_val(i+1, h+1, subQ(i+1,h+1-j));
  			}
		  }
    }
    return Q;
}

Vector SD(Matrix& A, Vector& v, int& count){
	assert(nrows(A) == ncols(A));
	int m =nrows(A);
	// we want a square system

	// INITIALISATION
	Vector x= v;
	double tol = 1.e-14*v.norm();
	double top;
	Vector r = A*x;
	r= r-v;
	double bottom;
	Vector temp(m);
	double alpha;
	//std::cout<<" before the while loop "<<(A*x -c).norm() << std::endl<<std::flush;
	while (r.norm() > tol){
		top = r*r;
		temp = A*r;
		bottom = r*temp;
		alpha= top/bottom;
		count+=1;
		x = x + top*r/bottom;
		r = v - A*x;
	}
   //	std::cout<<count<<std::endl<<std::flush;

	return x;
}

Vector CG(Matrix& A, Vector& v, int& count){
	assert(nrows(A) == ncols(A));
	int m =nrows(A);
	// we want a square system

	// INITIALISATION
  Vector x= v;
	double tol = 1.e-14*v.norm();
	Vector r = v - A*x;
	Vector p = r;
	Vector temp(m);
	temp = A*r;
	double alpha;
	alpha = r*temp;
	alpha = ((r.norm())*(r.norm()))/alpha;
	x = x +alpha *r;
	r = v - A*x;
	double beta;
	double anorm;
	Vector temp2(m);
	//std::cout<<" before the while loop "<<(A*x -c).norm() << std::endl<<std::flush;
	while ((A*x -v).norm() > tol){
		temp = A*r;
		temp2 = A*p;
		beta = r*temp2;
		anorm = p*temp2;
		beta = beta /anorm;
		p = r-beta*p;
		alpha = r*p;
    temp2 = A*p;
    anorm = p*temp2;
		alpha = alpha/anorm;
		count+=1;
		x = x + alpha*p;
		r = v - A*x;
    if (count>1000000){
      count = -1;
      return x;
    }
	}
  //	std::cout<<count<<std::endl<<std::flush;
	return x;
}


Vector CG_pre(Matrix& B, Vector& b, int& count){
	assert(nrows(B) == ncols(B));
	int m =nrows(B);
	// we want a square system
  // do the preconditioning by multiplying with D^{-1}
  Matrix A = B;
  for (int k=1;k<=m;k++){
    double x = A(k,k);
    //A.set_col(k, A.get_col(k)/x);
    for (int j=1;j<=m;j++){
     A.set_val(k,j, A(k,j)/x);
    }
  }
  Vector v = b;
  for (int k=1;k<=m;k++){
    v.set_val(k, v(k)/B(k,k));
  }
	// INITIALISATION
	Vector x= v;
	double tol = 1.e-12;
	Vector r = v - A*x;
	Vector p = r;
	Vector temp(m);
	temp = A*r;
	double alpha;
	alpha = r*temp;
	alpha = ((r.norm())*(r.norm()))/alpha;
	x = x +alpha *r;
	r = v - A*x;
	double beta;
	double anorm;
	Vector temp2(m);
	//std::cout<<" before the while loop "<<(A*x -c).norm() << std::endl<<std::flush;
	while (r.norm() > tol){
    count+=1;
		temp = A*r;
		temp2 = A*p;
		beta = r*temp2;
		anorm = p*temp2;
		beta = beta /anorm;
		p = r-beta*p;
		alpha = r*p;
    temp2 = A*p;
    anorm = p*temp2;
		alpha = alpha/anorm;

		x = x + alpha*p;
		r = v - A*x;
    if (count>100){
      r = r*0;
      count = -100000;
    }
  }
  //	std::cout<<count<<std::endl<<std::flush;
	return x;
}

Vector momentum(Matrix& A, Vector& v, int& count,double lambda,double mu){
	assert(nrows(A) == ncols(A));
	int m =nrows(A);
	// we want a square system

	// INITIALISATION
	Vector x= v;
	double tol = 1.e-12;
	Vector r = A*x;
	r= r-v;
	Vector temp(m);
	//std::cout<<" before the while loop "<<(A*x -c).norm() << std::endl<<std::flush;
	while ((A*x -v).norm() > tol){
		count+=1;
    temp = mu*temp - lambda*r;
		x = x + temp;
		r = v - A*x;
	}
   //	std::cout<<count<<std::endl<<std::flush;

	return x;
}
