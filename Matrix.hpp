#ifndef MATRIXDEF
#define MATRIXDEF



//  **********************
//  *  Class of Matrix  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written


#include <cmath>
#include "Exception.hpp"//  This class throws errors using the class "error"
#include "Vector.hpp"

class Matrix
{
private:
   // member variables
   double** mData;   // data stored in vector of vectors
   int mShape[2];      // shape of matrix

public:

   Matrix(int rows, int cols);
   Matrix(const Matrix& C);
   double& operator()(int i, int j);
   Matrix& operator=(const Matrix& v);
   Matrix T() const; // transpose
   // Unary operator
   friend Matrix operator-(const Matrix& v);

   // destructor
   ~Matrix();
   //getting shape of matrix
   friend int nrows(const Matrix& C);
   friend int ncols(const Matrix& C);
   //get col or row vector as double*
   Vector get_col(int k);
   Vector get_row(int k);
   void set_val(int i, int j, double x);
   //void set_col(int i, Vector& x);

   // add and substract
   friend Matrix operator+(const Matrix& v1, const Matrix& v2);
   friend Matrix operator+(const Matrix& A, const double k);
   friend Matrix operator-(const Matrix& v1, const Matrix& v2);
   // Multiplications
   friend Matrix operator*(Matrix& v1, Matrix& v2);
   friend Vector operator*(Matrix& A, Vector& b);
   friend Matrix operator*(const Matrix& v, const double& a);
   friend Matrix operator*(const double& a, const Matrix& v);
   friend Matrix operator/(const double& a, const Matrix& v);
   Vector diag();
   friend Matrix diag(const Matrix& C);
   friend std::ostream& operator<<(std::ostream& output, const Matrix& C);
};
int nrows(const Matrix& C);
int ncols(const Matrix& C);
Matrix operator-(const Matrix& v);
Matrix diag(const Matrix& C);
Matrix diag( Vector& v);
Matrix operator+(const Matrix& v1, const Matrix& v2);
Matrix operator+(const Matrix& A, const double k);
Matrix operator-(const Matrix& v1, const Matrix& v2);

//double s_prod(Vector& a, Vector& b,int length);
//Matrix operator*(const Matrix& v1, const Matrix& v2);
Matrix operator*( Matrix& v1,  Matrix& v2);
Matrix operator*(const Matrix& v, const double& a);
Matrix operator*(const double& a, const Matrix& v);
Matrix operator/(const double& a, const Matrix& v);
Vector operator*(Matrix& A, Vector& b);

Vector Jacobi(Matrix& A, Vector& v, int& count);
Vector SOR(Matrix& A, Vector& v, double w, int& count);
Vector SD(Matrix& A, Vector& v, int& count);
Vector momentum(Matrix& A, Vector& v, int& count, double lambda,double mu);

Vector CG(Matrix& A, Vector& v, int& count);
Vector CG_pre(Matrix& A, Vector& v, int& count);

Vector LUsolve(Matrix& A, Vector& v);
Vector QRsolve(Matrix& A, Vector& v);
Vector QRsolve2( Matrix& A,  Vector& v);
Matrix randm(Matrix& A);
Matrix qr_q(Matrix& A);
Matrix rand_basis_gs(Matrix& A);
Matrix rand_basis_angle(Matrix& A);


// graph functions
Matrix er_graph(Matrix& A,float p);
Matrix floyd_war(Matrix& A);
Vector evec_c(Matrix& A);

#endif
