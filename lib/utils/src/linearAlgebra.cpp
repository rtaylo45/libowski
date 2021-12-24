//*****************************************************************************
// Author: Zack Taylor
//
// Linear algebra utilities that Eigen doesn't provide
//*****************************************************************************
#include "linearAlgebra.h"
#include <iostream>

//*****************************************************************************
// Implements the Moore-Penrose pseudo inverse for a eigen sparse matrix
//
// @param A    matrix
//*****************************************************************************
template <typename derived>
SparseMatrix<derived> MoorePenroseInv(const SparseMatrix<derived>& A){
  SparseMatrix<derived> Ainv;
  SparseMatrix<derived> AConjugateTranspose;
  SparseMatrix<derived> Atemp;
  SparseMatrix<derived> I(A.rows(),A.cols());
  SparseQR<SparseMatrix<derived>, COLAMDOrdering<int>> solver;

  I.setIdentity();
  AConjugateTranspose = A.adjoint();
  Atemp = AConjugateTranspose*A;
  solver.compute(Atemp);
  Ainv = solver.solve(I)*AConjugateTranspose;
  return Ainv;
}

//*****************************************************************************
// l1 norm of a sparse matrix
//
// @param A    Sparse matrix
//*****************************************************************************
double l1norm(const SparseMatrixD& A){
  double norm;
  norm = (A.cwiseAbs()*VectorD::Ones(A.cols())).maxCoeff();
  return norm;

}

//*****************************************************************************
// Produces the l1norm of A*B. Need to eventually change this to use the
// estiment of the norm which is use in the paper.
//
// @param A    Sparse matrix
// @param B    Sparse matrix
//*****************************************************************************
double normAm(const SparseMatrixD& A, const SparseMatrixD& B){
  SparseMatrixD C = A*B;
  double C1norm = l1norm(C);
  return C1norm;
}

//*****************************************************************************
// Produces the 1lnorm of A^m
//
// @param A    Sparse matrix
// @param m    inteter, power of the matrix
//*****************************************************************************
double normAm(const SparseMatrixD& A, const int m){
  SparseMatrixD C = A;
  double C1norm;

  // Rises the matrix to the power m
  for (int i; i < m; i++){
    C = C*A;
  }
  C1norm = l1norm(C);
  return C1norm;
}

//*****************************************************************************
// Arnoldi algorithm
//
// Used for generating the Krylov subspace. This process produces and
// orthnormal basis V(m x n) and an upper Hessenberg matrix (n x n). The
// algorithm produces coefficinet for h(n+1,n) and vector v(n+1) but these
// are not added to the returned matrices. In the future these can be added
// to produced to more accurate Krylov subspace approximation.
//
// The orthnormal basis is a dense matrix and will be retured as one H is
// lower triangular so it is retured as sparse. H is also the matrix that
// will be evaluated with the matrix exponential so it needs to be sparse to
// be evaluated by the existing matrix exp methods.
//
// @param A    Sparse matrix
// @param b    b vector used to build the subspace
// @param n    Dimention of the subspace
// @param V    Returned orthonormal basis
// @param H    Returned upper Hessenburg matrix, the projection of A into the
//          Krylov subspace
//*****************************************************************************
void arnoldi(const SparseMatrixD& A, const VectorD& b, const int n,
  MatrixD& V, SparseMatrixD& H){
  VectorD q, v;
  double h, normV;
  const int cols = A.cols();
  const int rows = A.rows();
  const double eps = 1.e-12;

  // if the subspace dim is greater than the A matrix dim, kill program
  if (n > cols) {
    std::string errorMessage = " The Krylov subspace dimension is \n"
      " greater than the input matrix size\n";
    libowskiException::runtimeError(errorMessage);
  }

  // Set temp matrices to build the subspace
  MatrixD Vtemp = MatrixD::Zero(cols, n+1);
  SparseMatrixD Htemp(n+1, n);
  // reserve space for the H matrix. This includes the main diagonal
  // the sub main diagonal, and all elements above the main diagonal.
  Htemp.reserve(n*(n+1) - (n-1)*n/2 - 1);

  q = b/b.norm();
  Vtemp.col(0) = q;

  // Loop over columns
  for (int i = 0; i < n; i++){
    // new candidant vector
    v = A*q;
    // Subtract previous vectors
    for (int j = 0; j <= i; j++){
      h = Vtemp.col(j).conjugate().dot(v);
      Htemp.insert(j, i) = h;
      v = v - h*Vtemp.col(j);

    }
    normV = v.norm();
    Htemp.insert(i+1,i) = normV;
    // If v is too small its the zero vector and should not be added
    if (normV > eps){
      q = v/normV;
      Vtemp.col(i+1) = q;
    }
    else{
      V = Vtemp.topLeftCorner(rows,n);
      H = Htemp.topLeftCorner(n,n);
      return;
    }
  }
  V = Vtemp.topLeftCorner(rows,n);
  H = Htemp.topLeftCorner(n,n);
  return;
}

// Data types that can use the template functions
template SparseMatrixLD MoorePenroseInv(const SparseMatrixLD& A);
template SparseMatrixD MoorePenroseInv(const SparseMatrixD& A);
template SparseMatrixCLD MoorePenroseInv(const SparseMatrixCLD& A);
template SparseMatrixCD MoorePenroseInv(const SparseMatrixCD& A);
