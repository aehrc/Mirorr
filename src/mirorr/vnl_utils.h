/*=========================================================================
Program: mirorr
Module: vnl_utils.h
Author: David Rivest-Henault
Created: 26 Oct 2012

Copyright (c) 2009-15 CSIRO. All rights reserved.

For complete copyright, license and disclaimer of warranty
information see the LICENSE.txt file for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef VNL_UTILS_H_
#define VNL_UTILS_H_

#ifdef DEBUG_VNL_UTILS
#include <iostream>
#endif

#include <limits>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_matrix_inverse.h>

#if ITK_VERSION_MAJOR < 4
#include <vnl/algo/vnl_svd.h>
#else
#include <vnl/algo/vnl_svd_fixed.h>
#endif



/**
 * Given a positive definite matrix M, compute M^(1/2),
 * such that M = M^(1/2) * M^(1/2) using Denman–Beavers iterations.
 *
 * @see http://en.wikipedia.org/wiki/Matrix_square_root#By_Denman.E2.80.93Beavers_iteration
 * @see Reuter, Rosas, Fischl. "Highly accurate inverse consistent registration:
 * A robust approach", NeuroImage (53), 2010. Equation (29), p.1188.
 */

template<class TVnlMatrix>
TVnlMatrix
vnl_matrix_sqrt(const TVnlMatrix &M, const double EPS = 1e-10,
      const size_t MAX_ITER=20)
{

  if (M.cols() != M.rows()) {
    return TVnlMatrix(M);
  }

  typedef typename TVnlMatrix::element_type T;
  vnl_matrix<T> Y = M;
  vnl_matrix<T> Z(M.rows(),M.cols());
  Z.set_identity();

  // iterative delta
  vnl_matrix<T> D;
  double d = std::numeric_limits<double>::max();

  size_t i = 0;
  vnl_matrix<T> newY, newZ;
  while (d > EPS && i < MAX_ITER)
  {
    newY = 0.5 * (Y + vnl_matrix_inverse<T>(Z));
    newZ = 0.5 * (Z + vnl_matrix_inverse<T>(Y));

    Y = newY;
    Z = newZ;

    i++;
    D = (Y * Y) - M;
    d = D.absolute_value_max();

  }

  if (i >= MAX_ITER && d > EPS*1e3) { //Only display if the error is really big
      std::cout << "WARNING: reached the maximum number "
          << "of iteration without achieving the requested precision of "
          << EPS << "\n"
          << "         In vnl_matrix_sqrt(...), current precision is "
          << d << std::endl;
  }

#ifdef DEBUG_VNL_UTILS
  std::cout << "Computed vnl_matrix_sqrt in " << i << " iterations"
      << " with a minimal precision of " << d << std::endl;
#endif

  return Y;
}

/**
 * Compute the square root of the matrix m12 = sqrt(M)
 * and im12 the inverse of m12 using im12 = m12 * inv(M) for consistency
 */
template<class TVnlMatrix>
bool
vnl_matrix_sqrt_and_inverse_sqrt(const TVnlMatrix &M,
    TVnlMatrix &m12, TVnlMatrix &im12,
    const double EPS = 1e-10, const size_t MAX_ITER=20)
{
  typedef typename TVnlMatrix::element_type T;

  if (M.cols() != M.rows()) {
    return false;
  }

  m12 = vnl_matrix_sqrt(M,EPS,MAX_ITER);

  //im12 = m12 * vnl_matrix_inverse<T>(M); //Preserve overall transform better
  im12 = vnl_matrix_inverse<T>(m12); //Preserve symmetry better

  return true;
}

/**
 * Compute an orthogonal matrix that is 'close' to the input matrix
 */
template<class T, int N>
vnl_matrix_fixed<T,N,N>
orthogonalize(vnl_matrix_fixed<T,N,N> m)
{
#if ITK_VERSION_MAJOR < 4
  vnl_svd<T> svd(m);
#else
  vnl_svd_fixed<T,N,N> svd(m);
#endif
  
  return svd.U() * svd.V().conjugate_transpose();
}


/**
 * Copy the 3x3 top left sub-matrix of H in M
 */
template<class TVnlMatrixA, class TVnlMatrixB>
bool
copy3x3SubMatrix(const TVnlMatrixA &H, TVnlMatrixB *M)
{
  typedef typename TVnlMatrixB::element_type T;

  if (H.cols() < 3 || H.rows() < 3 || M->cols() < 3 || M->rows() < 3) {
    return false;
  }

  TVnlMatrixB &m = *M;
  m(0,0) = static_cast<T>(H(0,0)); m(0,1) = static_cast<T>(H(0,1)); m(0,2) = static_cast<T>(H(0,2));
  m(1,0) = static_cast<T>(H(1,0)); m(1,1) = static_cast<T>(H(1,1)); m(1,2) = static_cast<T>(H(1,2));
  m(2,0) = static_cast<T>(H(2,0)); m(2,1) = static_cast<T>(H(2,1)); m(2,2) = static_cast<T>(H(2,2));

  return true;
}


#endif /* VNL_UTILS_H_ */


