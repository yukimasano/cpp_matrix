# cpp_matrix: A C++ implementation of a matrix class and linear solvers 
(For code description see bottom)

This project's purpose was to analyze the different scaling patterns of commonly used Linear solvers.
The algorithms implemented here are: LU, QR (implicit and explicit), SD, CG, CG with pre-conditioning (pc), Sucessive-overrelaxation, Jacobi, Gauss-Seidel.

# Linear Systems?

Linear systems have the shape *Ax=b* with a known matrix A and known vector b. Essentially all we're trying to do is `x=A\b`, in Matlab notation. Note also that these solvers can be used for finding an x that minimizes *| Ax - b |*. This is done by solving the symmetric system *A'Ax = A'b* (note however that the condition number will be squared).

# Quick overview of when which solver can be used.

## Any A

LU, QR, SD, (if you take the squared system: CG)

## Symmetric A

LU, QR, SD, CG.

## Strictly-row/column-diagonally-dominant A

LU, QR, SD, CG, Sucessive-overrelaxation (includes Jacobi and Gauss-Seidel)

# Results

## Same condition number κ, varying size N, non-symmetric A

## Same condition number κ, varying size N, symmetric A

## Same size number N, condition number κ, symmetric A

## SRDD with varying size N = κ, symmetric A
