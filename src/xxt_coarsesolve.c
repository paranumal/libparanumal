#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "crs.h"
#include "mpi.h"

void * xxt_setup(uint num_local_rows, 
                 uint* row_ids,
                 uint nnz, 
                 uint*   A_i,
                 uint*   A_j,
                 double* A_vals,
                 int null_space) {
  int np, my_id;
  struct comm COMM;
  struct crs_data *crs_A;

  MPI_Comm_size(MPI_COMM_WORLD,&np);
  comm_init(&COMM,(comm_ext) MPI_COMM_WORLD);

  my_id = COMM.id;

  crs_A = crs_setup(num_local_rows, row_ids,
                  nnz, A_i, A_j, A_vals,
                  null_space, &COMM);
  crs_stats(crs_A);

  return (void *) crs_A;
}

int xxt_solve(double* x,
              void* crs_A,
              double* rhs) {
  struct crs_data *A;

  A = (struct crs_data*) crs_A;

  crs_solve(x,A,rhs);

  return 0;
}

int xxt_free(void* crs_A) {
  struct crs_data* A;
  A = (struct crs_data*) crs_A;

  crs_free(A);

  return 0;
}
