# cpp_matrix: A C++ implementation of a matrix class and linear solvers 
(For code description see bottom)

This project's purpose was to analyze the different scaling patterns of commonly used Linear solvers.
The algorithms implemented here are: LU, QR (implicit and explicit), SD, CG, CG with Jacobi pre-conditioning, Sucessive-overrelaxation, Jacobi, Gauss-Seidel.

# Linear Systems?

Linear systems have the shape `Ax = b` with a known matrix A and known vector b. Essentially all we're trying to do is `x=A\b`, in Matlab notation. Note also that these solvers can be used for finding an x that minimizes `|| Ax - b ||`. This is done by solving the symmetric system `A'Ax = A'b` (note however that the condition number will be squared).

# Quick overview of when which solver can be used.

## Any A

LU, QR, SD, (if you take the squared system: CG)

## Symmetric A

LU, QR, SD, CG.

## Strictly-row/column-diagonally-dominant A

LU, QR, SD, CG, Sucessive-overrelaxation (includes Jacobi and Gauss-Seidel)

# Results
Steps for generating random matrices of size `N` with desired condition numbers `κ`: 
 1. We generate random matrix Q of size `NxN`
 2. We orthonormalize the columns via modified Gram-Schmidt
 3. We let D = diag(linspace(1,...,`κ`))
 4. `A = QD` or if we want a symmetric system: `A = QDQ'`
 5. Generate random vector `x` and let `b = Ax`.
Steps for generating SRDD matrices with  size `N` and approximately `κ = 2N`
 1. Generate random matrix A of size `NxN`
 2. Set diagonal entries of A to linspace(N,...,2N)
 3. Generate random vector `x` and let `b = Ax`.
By [Gershogorin discs](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem) the matrix has approximately the eigenvalues (N,...,2N) and thus `κ = 2`.

## Same condition number `κ`, varying size `N`, non-symmetric A

<img src="https://user-images.githubusercontent.com/29401818/33230601-6e5347b2-d1de-11e7-93c8-a0e591901fda.png" height ="200">
In the above figure we use the normal system for the CG solvers.

## Same condition number `κ`, varying size `N`, symmetric A

<img src="https://user-images.githubusercontent.com/29401818/33230603-6e914bf2-d1de-11e7-9e83-02dd4116353b.png" height ="200">

## Same size `N`, varying condition number `κ`, symmetric A

<img src="https://user-images.githubusercontent.com/29401818/33230600-6e40a4a4-d1de-11e7-8589-9b668874fb42.png" height ="200">

## SRDD with condition number `κ=2` and varying size `N`

<img src="" height ="200">

# Critical remarks
 * Definitely easier to implement in Matlab, Python..
 * Matrices were too easy for my solvers in some cases by construction (CG with pre-conditioning converges in 1-step for `A = QDQ'`)
 * Writing everything in C++ was instructive but one should really use proper templates like [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
 
# Code
 * `vary_N.cpp`: implements the above procedure for symmetric A system with 10 iterations at every size. Sizes grow exponentially.
 * `vary_kappa.cpp`: same as vary_N except that now we vary `κ` at every iteration and keep `N = 100`
 * `SRDD.cpp`: Generates SRDD matrix as above and solves it, also 10 iterations at every size that grows exponentially. 
 * `Matrix.cpp`: Includes the matrix class and all methods and the solvers
 * `Vector.cpp`: A vector class, as from introductory C++ course of [Joe Pit-Francis](https://www.cs.ox.ac.uk/people/joe.pitt-francis/)
