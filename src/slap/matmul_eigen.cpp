//
// Created by Sam Schoedel on 5/2/23.
// Copyright (c) 2022 Robotic Exploration Lab. All rights reserved.
//

#include <iostream>

#include "Eigen/Dense"

// using namespace Eigen;

extern "C" {
#include "matmul.h"
#include "tri.h"
#include "copy_matrix.h"

enum slap_ErrorCode slap_MatMulAdd(
    Matrix C, Matrix A, Matrix B, sfloat alpha,
    sfloat beta) {
  
  // std::cout << "this is the eigen implementation" << std::endl;

  // // Check for special structure
  // if (slap_GetType(A) == slap_TRIANGULAR_UPPER) {
  //   return slap_UpperTriMulAdd(C, A, B, alpha, beta);
  // }
  // if (slap_GetType(A) == slap_TRIANGULAR_LOWER) {
  //   return slap_LowerTriMulAdd(C, A, B, alpha, beta);
  // }

  // SLAP_ASSERT_VALID(C, SLAP_INVALID_MATRIX, "MatMulAdd: invalid C matrix");
  // SLAP_ASSERT_VALID(A, SLAP_INVALID_MATRIX, "MatMulAdd: invalid A matrix");
  // SLAP_ASSERT_VALID(B, SLAP_INVALID_MATRIX, "MatMulAdd: invalid B matrix");

  int n = slap_NumRows(A);
  int m = slap_NumCols(A);
  int p = slap_NumCols(B);
  // SLAP_ASSERT(slap_NumRows(B) == m, SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             "MatMulAdd: dimension mismatch, B has %d rows, expected %d", slap_NumRows(B),
  //             m);
  // SLAP_ASSERT(slap_NumRows(C) == n, SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             "MatMulAdd: dimension mismatch, C has %d rows, expected %d", slap_NumRows(C),
  //             n);
  // SLAP_ASSERT(slap_NumCols(C) == p, SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             SLAP_INCOMPATIBLE_MATRIX_DIMENSIONS,
  //             "MatMulAdd: dimension mismatch, C has %d columns, expected %d",
  //             slap_NumCols(C), p);
              
  // MatrixXd A_eigen = Map<Matrix<sfloat, Dynamic, Dynamic>>(A.data, m, n);
  // MatrixXd B_eigen = Map<Matrix<sfloat, Dynamic, Dynamic>>(B.data, n, p);
  // MatrixXd C_eigen = Map<Matrix<sfloat, Dynamic, Dynamic>>(C.data, m, p);

  // C_eigen = beta * C_eigen + alpha * A_eigen * B_eigen;

  // slap_CopyFromArray(C, C_eigen.data());

  Eigen::Map<Eigen::Matrix<sfloat, Eigen::Dynamic, Eigen::Dynamic>>(C.data, m, p) = 
      beta * Eigen::Map<Eigen::Matrix<sfloat, Eigen::Dynamic, Eigen::Dynamic>>(C.data, m, p)
      + alpha * Eigen::Map<Eigen::Matrix<sfloat, Eigen::Dynamic, Eigen::Dynamic>>(A.data, m, n) 
      * Eigen::Map<Eigen::Matrix<sfloat, Eigen::Dynamic, Eigen::Dynamic>>(B.data, n, p);


  // for (int i = 0; i < n; ++i) {    // rows of output
  //   for (int j = 0; j < p; ++j) {  // Columns of output
  //     sfloat* Cij = slap_GetElement(C, i, j);
  //     *Cij *= beta;
  //     for (int k = 0; k < m; ++k) {  // columns of A, rows of B
  //       sfloat Aik = *slap_GetElementConst(A, i, k);
  //       sfloat Bkj = *slap_GetElementConst(B, k, j);
  //       *Cij += alpha * Aik * Bkj;
  //     }
  //   }
  // }
  return SLAP_NO_ERROR;
}

} // extern "C"
