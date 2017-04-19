#include "parAlmond.h"


void axpy(almond_t *almond, ell *A, dfloat alpha, occa::memory o_x, dfloat beta, occa::memory o_y) {

  if(A->Nrows){
    const iint numBlocks = (A->Nrows+AGMGBDIM-1)/AGMGBDIM;

    almond->ellAXPYKernel(numBlocks, AGMGBDIM, A->Nrows, A->nnzPerRow, A->strideLength, alpha, beta,
	                 A->o_cols, A->o_coefs, o_x, o_y);
  }
}

void zeqaxpy(almond_t *almond, ell *A, dfloat alpha, occa::memory o_x, 
            dfloat beta, occa::memory o_y,  occa::memory o_z) {

  if(A->Nrows){
    const iint numBlocks = (A->Nrows+AGMGBDIM-1)/AGMGBDIM;

    almond->ellZeqAXPYKernel(numBlocks, AGMGBDIM, A->Nrows, A->nnzPerRow, A->strideLength, alpha, beta,
	                   A->o_cols, A->o_coefs, o_x, o_y, o_z);
  }
}

