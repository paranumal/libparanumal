/* compile with C compiler (not C++) */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define UNDERSCORE 1
#define USE_NAIVE_BLAS 
#define NO_NEX_EXITT 1
#define GLOBAL_LONG_LONG 1
#define PREFIX jl_

#define MPI 1

#include "gslib.h"

void meshParallelGatherScatter(void *gsh, float *v){

  /* need gs_float or gs_double */
  gs(v, gs_float, gs_add, 0, gsh, 0);  

}


