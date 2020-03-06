#include "blas2.h"
#include "mkl_types.h"
#include "mkl_cblas.h"

void SGEMV(char TRANS,int M,int N,float ALPHA,float* A,int LDA,float* X,int INCX,float BETA,float* Y,int INCY)
{
      CBLAS_ORDER  order=CblasRowMajor;
      CBLAS_TRANSPOSE trans;
      if(TRANS=='n'||TRANS=='N')
        trans= CblasNoTrans;
      else
        trans= CblasTrans;

      cblas_sgemv(order, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    return;
}



