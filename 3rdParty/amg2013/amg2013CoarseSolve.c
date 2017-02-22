#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_sstruct_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "sstruct_mv.h"

typedef struct {
  hypre_ParCSRMatrix    *A;
  hypre_ParVector       *x;
  hypre_ParVector       *rhs;
  HYPRE_Solver          solver;
  HYPRE_Solver          precond;

  HYPRE_BigInt global_size;
  HYPRE_BigInt * row_starts;

  int numLocalRows;

  int     diag_nnz;
  int    *diag_i;      //local crs sparse matrix (locally indexed)
  int    *diag_j;
  double *diag_data;

  int     offd_nnz;
  int    *offd_i;      //nonlocal crs sparse matrix (globally indexed)
  int    *offd_j;     
  double *offd_data;
  HYPRE_BigInt *colMap;

  double *xDouble;
  double *rhsDouble;

  char* iintType;
  char* dfloatType;

} amg2013_t;


void * amg2013SetupCSR(int *row_starts,     //[numproc+1] global partition array   
                       int    *diag_i,      //local crs sparse matrix (locally indexed)
                       int    *diag_j,
                       void   *diag_data,
                       int    *offd_i,      //nonlocal crs sparse matrix (globally indexed)
                       int    *colMap,
                       void   *offd_data,
                       const char* iintType, 
                       const char* dfloatType) {
  //defaults
  double tol = 1.e-6;
  int    maxit_prec = 100;
  int    maxit_sol = 500;
  
  int num_procs, myid;
  int numLocalRows;
  int n,i,j;

  hypre_ParCSRMatrix *A;
  hypre_CSRMatrix *diag;
  hypre_CSRMatrix *offd;

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  amg2013_t *amg = (amg2013_t*) calloc(1, sizeof(amg2013_t));


  if (!strcmp(iintType,"long")) { 
    printf("Exiting due to use of ulong in iintType %s\n", iintType);
    exit(-1);
  }

  if (!strcmp(dfloatType,"float")) { 
    amg->xDouble   = (double *) malloc(numLocalRows*sizeof(double));
    amg->rhsDouble = (double *) malloc(numLocalRows*sizeof(double));
  } 

  amg->global_size = (HYPRE_BigInt) row_starts[num_procs];
  amg->row_starts = (HYPRE_BigInt*) malloc((num_procs+1)*sizeof(HYPRE_BigInt));
  for (n =0;n<=num_procs;n++) amg->row_starts[n] = (HYPRE_BigInt) row_starts[n];
  
  numLocalRows = row_starts[myid+1]-row_starts[myid];
  amg->numLocalRows = numLocalRows;
  amg->diag_nnz = diag_i[numLocalRows];
  amg->offd_nnz = offd_i[numLocalRows];

  amg->diag_i = (int *) malloc((numLocalRows+1)*sizeof(int));
  amg->diag_j = (int *) malloc(amg->diag_nnz*sizeof(int));
  amg->diag_data = (double *) malloc(amg->diag_nnz*sizeof(double));

  amg->offd_i = (int *) malloc((numLocalRows+1)*sizeof(int));
  if (amg->offd_nnz) {
    amg->offd_j = (int *) malloc(amg->offd_nnz*sizeof(int));
    amg->colMap = (HYPRE_BigInt *) malloc(amg->offd_nnz*sizeof(HYPRE_BigInt));
    amg->offd_data = (double *) malloc(amg->offd_nnz*sizeof(double));
  } else {
    amg->offd_j = NULL;
    amg->colMap = NULL;
    amg->offd_data = NULL;
  }

  //copy data into amg struct
  for (n=0;n<numLocalRows+1;n++) {
    amg->diag_i[n] = diag_i[n];
    amg->offd_i[n] = offd_i[n];
  }
  for (n=0;n<amg->diag_nnz;n++) {
    amg->diag_j[n] = diag_j[n];
    if (!strcmp(dfloatType,"float")) 
      amg->diag_data[n] = ((float*) diag_data)[n];
    else 
      amg->diag_data[n] = ((double*) diag_data)[n];
  }
  for (n=0;n<amg->offd_nnz;n++) {
    amg->colMap[n] = (HYPRE_BigInt) colMap[n];
    if (!strcmp(dfloatType,"float")) 
      amg->offd_data[n] = ((float*) offd_data)[n];
    else 
      amg->offd_data[n] = ((double*) offd_data)[n];
  }


  //allocate 
  A = hypre_ParCSRMatrixCreate(MPI_COMM_WORLD, 
                               amg->global_size, amg->global_size, //Global matrix size
                               amg->row_starts, amg->row_starts,   //numproc+1 size global partition array
                               amg->offd_nnz,   //total number of nonlocal nonzero colunms 
                               amg->diag_nnz,   //total number of local nonzero entries   
                               amg->offd_nnz);  //total number of nonlocal nonzero entries 

  amg->rhs = hypre_ParVectorCreate(MPI_COMM_WORLD, amg->global_size, amg->row_starts);
  amg->x = hypre_ParVectorCreate(MPI_COMM_WORLD, amg->global_size, amg->row_starts);

  amg->A = A;

  diag = hypre_ParCSRMatrixDiag(A);
  hypre_CSRMatrixI(diag) = amg->diag_i;
  hypre_CSRMatrixJ(diag) = amg->diag_j;
  hypre_CSRMatrixData(diag) = amg->diag_data;

  int * tmp_j;
  if (amg->offd_nnz) {
    tmp_j = (int*) malloc(amg->offd_nnz*sizeof(int));  
  }
  for (i=0; i < amg->offd_nnz; i++) {
    amg->offd_j[i] = i;
    tmp_j[i] = i;
  }
  if (num_procs > 1) {
    hypre_BigQsortbi(amg->colMap, tmp_j, 0, amg->offd_nnz-1);
    for (i=0; i < amg->offd_nnz; i++){
      for (j=0; j < amg->offd_nnz; j++){
        if (amg->offd_j[i] == tmp_j[j]){
          amg->offd_j[i] = j;
          break;
        }
      }
    }
  }
  

  offd = hypre_ParCSRMatrixOffd(A);
  hypre_CSRMatrixI(offd) = amg->offd_i;
  if (amg->offd_nnz) {
    hypre_CSRMatrixJ(offd) = amg->offd_j;
    hypre_CSRMatrixData(offd) = amg->offd_data;
    hypre_ParCSRMatrixColMapOffd(A) = amg->colMap;
    free(tmp_j);   
  }
  
  amg->iintType = strdup(iintType);
  amg->dfloatType = strdup(dfloatType);

  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &(amg->solver));
  HYPRE_PCGSetTol( amg->solver, tol );
  HYPRE_PCGSetTwoNorm( amg->solver, 1 );
  HYPRE_PCGSetRelChange( amg->solver, 0 );
  HYPRE_PCGSetPrintLevel( amg->solver, 0 );
  

  /* use BoomerAMG as preconditioner */
  HYPRE_PCGSetMaxIter( amg->solver, maxit_prec);
  HYPRE_BoomerAMGCreate(&(amg->precond)); 
  HYPRE_BoomerAMGSetCoarsenType(amg->precond, 10);
  HYPRE_BoomerAMGSetStrongThreshold(amg->precond, 0.25);
  HYPRE_BoomerAMGSetAggNumLevels(amg->precond, 1);
  HYPRE_BoomerAMGSetInterpType(amg->precond, 6);
  HYPRE_BoomerAMGSetPMaxElmts(amg->precond, 4);
  HYPRE_BoomerAMGSetTol(amg->precond, 0.0);
  HYPRE_BoomerAMGSetRelaxType(amg->precond, 8);
  HYPRE_BoomerAMGSetCycleRelaxType(amg->precond, 8, 3);
  HYPRE_BoomerAMGSetCycleNumSweeps(amg->precond, 1, 3);
  HYPRE_BoomerAMGSetRelaxOrder(amg->precond, 0);
  HYPRE_BoomerAMGSetPrintLevel(amg->precond, 0);
  HYPRE_BoomerAMGSetPrintFileName(amg->precond, "sstruct.out.log");
  HYPRE_BoomerAMGSetMaxIter(amg->precond, 1);
  HYPRE_PCGSetPrecond( amg->solver,
                (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                amg->precond );
  HYPRE_PCGSetup( amg->solver, (HYPRE_Matrix) amg->A,
             (HYPRE_Vector) amg->rhs, (HYPRE_Vector) amg->x );

  return (void*) amg;
}


void * amg2013SetupCOO(int *row_starts,     //[numproc+1] global partition array 
                    int    diag_nnz, 
                    int    *Ai,      //local coo sparse matrix (locally indexed)
                    int    *diag_j,
                    void   *diag_data,
                    int    offd_nnz,
                    int    *Bi,      //nonlocal coo sparse matrix (globally indexed)
                    int    *colMap,
                    void   *offd_data,
                    const char* iintType, 
                    const char* dfloatType) {

  void *amg;

  int num_procs, myid;
  int numLocalRows;
  int *diag_i, *offd_i;
  int n;


  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  numLocalRows = row_starts[myid+1]-row_starts[myid];

  diag_i = (int*) calloc(numLocalRows+1,sizeof(int));
  offd_i = (int*) calloc(numLocalRows+1,sizeof(int));


  //Contruct CSR i vector
  diag_i[0] = 0;
  offd_i[0] = 0;
  int cnt   = 0;
  int o_cnt = 0;

  int rowIndex =1;
  for (n=0;n<diag_nnz;n++){
    if (Ai[n]>=rowIndex) {
      diag_i[rowIndex++] = cnt;
      n--;
    } else {
      cnt++;
    }
  }
  diag_i[numLocalRows] = diag_nnz;

  rowIndex =1;
  for (n=0;n<offd_nnz;n++) {
    if (Bi[n]>=rowIndex) {
      offd_i[rowIndex++] = o_cnt;
      n--;
    }
    else {
      o_cnt++;
    }
  }
  offd_i[numLocalRows] = offd_nnz;

  amg = amg2013SetupCSR(row_starts,    
                         diag_i, diag_j, diag_data,
                         offd_i, colMap, offd_data,
                         iintType, dfloatType);

  free(diag_i);
  free(offd_i);

  return amg; 
}

int amg2013Solve(void* x,
                 void* AMG,
                 void* rhs) {
  
  int    num_iterations;
  double final_res_norm;

  int n;

  hypre_Vector *rhsLocal;
  hypre_Vector *xLocal;

  amg2013_t *amg = (amg2013_t *) AMG;

  if (!strcmp(amg->dfloatType,"float")) { 
    for (n=0;n<amg->numLocalRows;n++) {
      amg->xDouble[n]   = ((float*) x)[n];
      amg->rhsDouble[n] = ((float*) rhs)[n];
    }
  } else { //double
    amg->xDouble   = (double*) x;
    amg->rhsDouble = (double*) rhs;
  }

  xLocal   = hypre_ParVectorLocalVector(amg->x);
  rhsLocal = hypre_ParVectorLocalVector(amg->rhs);
  hypre_VectorData(xLocal)    = amg->xDouble;
  hypre_VectorData(rhsLocal) = amg->rhsDouble;


  HYPRE_PCGSolve(amg->solver, (HYPRE_Matrix) amg->A,
             (HYPRE_Vector) amg->rhs, (HYPRE_Vector) amg->x);
  
  
  HYPRE_PCGGetNumIterations( amg->solver, &num_iterations );
  HYPRE_PCGGetFinalRelativeResidualNorm( amg->solver, &final_res_norm );

  if (!strcmp(amg->dfloatType,"float")) { 
    for (n=0;n<amg->numLocalRows;n++) {
      ((float*) x)[n]   = (float) amg->xDouble[n];
      ((float*) rhs)[n] = (float) amg->rhsDouble[n];
    }
  }

  return 0;
}


int amg2013Free(void* AMG) {
  amg2013_t *amg = (amg2013_t *) AMG;

  HYPRE_ParCSRPCGDestroy(amg->solver);

  HYPRE_BoomerAMGDestroy(amg->precond);

  free(amg->diag_i);
  free(amg->diag_j);
  free(amg->diag_data);

  free(amg->offd_i);
  if (amg->offd_nnz) {
    free(amg->offd_j);
    free(amg->colMap);
    free(amg->offd_data);
  }

  hypre_ParCSRMatrixDestroy(amg->A);
  hypre_ParVectorDestroy(amg->rhs);
  hypre_ParVectorDestroy(amg->x); 

  if (!strcmp(amg->dfloatType,"float")) { 
    free(amg->xDouble);
    free(amg->rhsDouble);
  }

  return 0;
}