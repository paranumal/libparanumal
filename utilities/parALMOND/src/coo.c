#include "parAlmond.h"


void ax(almond_t *almond, coo *C, dfloat alpha, occa::memory o_x, occa::memory o_y) {

  const iint numBlocks = (C->nnz+AGMGBDIM-1)/AGMGBDIM;

  if(C->nnz){
    // do block-wise product
    almond->cooAXKernel1(numBlocks, AGMGBDIM, C->nnz, alpha, C->o_rows, C->o_cols, C->o_coefs,
	      o_x, o_y, C->o_temp_rows, C->o_temp_Ax);

    almond->cooAXKernel2(1, 1, numBlocks, C->o_temp_rows, C->o_temp_Ax, o_y);
  }
}

