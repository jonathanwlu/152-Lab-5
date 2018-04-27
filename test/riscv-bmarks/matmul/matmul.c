// See LICENSE for license details.

#include "dataset.h"
#include "util.h"
#include <stddef.h>

void matmul(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;

  for (i = 0; i < lda; i++) {
    for (j = coreid; j < lda; j += ncores) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
    }
  }
}

void matmul_opt(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;
  size_t a = lda / ncores;

  for (i = coreid * a; i < a + coreid * a; i = i + 2) {
    for (j = 0; j < lda; j = j + 2) {
      data_t sum00 = 0;
      data_t sum01 = 0;
      data_t sum10 = 0;
      data_t sum11 = 0;

      for (k = 0; k < lda; k++) {
        sum00 += B[k*lda + j] * A[i*lda + k];
        sum01 += B[k*lda + (j+1)] * A[i*lda + k];
        sum10 += B[k*lda + j] * A[(i+1)*lda + k];
        sum11 += B[k*lda + (j+1)] * A[(i+1)*lda + k];
      }

      C[j + i*lda] = sum00;
      C[(j+1) + i*lda] = sum01;
      C[j + (i+1)*lda] = sum10;
      C[(j+1) + (i+1)*lda] = sum11;
    }
  }
}