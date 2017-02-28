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
  HYPRE_ParCSRMatrix    A;
  HYPRE_ParVector       par_x;
  HYPRE_ParVector       par_rhs;
  HYPRE_Solver          solver;
  HYPRE_Solver          precond;
  
  hypre_ParVector       *x;
  hypre_ParVector       *rhs;

  HYPRE_BigInt global_size;
  HYPRE_BigInt *row_starts;

  int numLocalRows;
  int Nnum;
  int recvNnum;

  int     diag_nnz;
  int    *diag_i;      //local crs sparse matrix (locally indexed)
  int    *diag_j;
  double *diag_data;

  int     offd_nnz;
  int    *offd_i;      //nonlocal crs sparse matrix (globally indexed)
  int    *offd_j;     
  double *offd_data;
  HYPRE_BigInt *colMap;

  int *sendSortId, *globalSortId, *compressId;
  int *sendCounts, *sendOffsets,  *recvCounts, *recvOffsets;

  double *xAssembled;
  double *rhsAssembled;

  double *xUnassembled;
  double *rhsUnassembled;

  double *xSort;
  double *rhsSort;

  char* iintType;
  char* dfloatType;

} amg2013_t;


void * amg2013SetupCSR(int    Nnum,            //unassembled size (Nverts*Nelements)
                       int    *row_starts,     //[numproc+1] global partition array   
                       int    *diag_i,      //local crs sparse matrix (locally indexed)
                       int    *diag_j,
                       double *diag_data,
                       int    *offd_i,      //nonlocal crs sparse matrix (globally indexed)
                       int    *colMap,
                       double *offd_data,
                       int    *sendSortId, 
                       int    *globalSortId, 
                       int    *compressId,
                       int    *sendCounts, 
                       int    *sendOffsets, 
                       int    *recvCounts, 
                       int    *recvOffsets,
                       const char* iintType, 
                       const char* dfloatType) {
  //defaults
  double tol = 1.e-8;
  int    maxit_prec = 40000;
  int    maxit_sol = 500;
  
  int num_procs, myid;
  int numLocalRows;
  int n,i,j;

  hypre_ParCSRMatrix *A;
  hypre_Vector *rhsLocal;
  hypre_Vector *xLocal;
  
  hypre_CSRMatrix *diag;
  hypre_CSRMatrix *offd;

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  amg2013_t *amg = (amg2013_t*) calloc(1, sizeof(amg2013_t));

  amg->global_size = (HYPRE_BigInt) row_starts[num_procs];
  amg->Nnum = Nnum;
  amg->row_starts = (HYPRE_BigInt*) malloc((num_procs+1)*sizeof(HYPRE_BigInt));
  for (n =0;n<=num_procs;n++) amg->row_starts[n] = (HYPRE_BigInt) row_starts[n];

  numLocalRows = row_starts[myid+1]-row_starts[myid];
  amg->numLocalRows = numLocalRows;
  amg->recvNnum = compressId[numLocalRows];

  amg->sendSortId = (int*) calloc(Nnum,sizeof(int));
  amg->globalSortId = (int*) calloc(Nnum,sizeof(int));
  for (n=0;n<Nnum;n++) {
    amg->sendSortId[n] = sendSortId[n];
    amg->globalSortId[n] = globalSortId[n];
  }

  amg->sendCounts = (int*) calloc(num_procs,sizeof(int));
  amg->recvCounts = (int*) calloc(num_procs,sizeof(int));
  for (n=0;n<num_procs;n++){
    amg->sendCounts[n] = sendCounts[n]*sizeof(double);
    amg->recvCounts[n] = recvCounts[n]*sizeof(double);
  }

  amg->sendOffsets = (int*) calloc(num_procs+1,sizeof(int));
  amg->recvOffsets = (int*) calloc(num_procs+1,sizeof(int));
  for (n=0;n<num_procs+1;n++){
    amg->sendOffsets[n] = sendOffsets[n]*sizeof(double);
    amg->recvOffsets[n] = recvOffsets[n]*sizeof(double);
  }

  amg->compressId  = (int*) calloc(numLocalRows+1,sizeof(int));
  for (n=0;n<numLocalRows+1;n++) amg->compressId[n] = compressId[n];

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
    amg->diag_j[n]    = diag_j[n];
    amg->diag_data[n] = diag_data[n];
  }
  for (n=0;n<amg->offd_nnz;n++) {
    amg->colMap[n] = (HYPRE_BigInt) colMap[n];
    amg->offd_data[n] = offd_data[n];
  }

  amg->iintType = strdup(iintType);
  amg->dfloatType = strdup(dfloatType);

  int Nmax =  (Nnum >amg->recvNnum)? Nnum : amg->recvNnum;
  amg->xUnassembled   = (double *) malloc(Nmax*sizeof(double));
  amg->rhsUnassembled = (double *) malloc(Nmax*sizeof(double));
  amg->xSort   = (double *) malloc(Nmax*sizeof(double));
  amg->rhsSort = (double *) malloc(Nmax*sizeof(double));

  amg->xAssembled   = (double *) calloc(numLocalRows,sizeof(double));
  amg->rhsAssembled = (double *) calloc(numLocalRows,sizeof(double));


  //allocate 
  A = hypre_ParCSRMatrixCreate(MPI_COMM_WORLD, 
                               amg->global_size, amg->global_size, //Global matrix size
                               amg->row_starts, amg->row_starts,   //numproc+1 size global partition array
                               amg->offd_nnz,   //total number of nonlocal nonzero colunms 
                               amg->diag_nnz,   //total number of local nonzero entries   
                               amg->offd_nnz);  //total number of nonlocal nonzero entries 

  amg->rhs = hypre_ParVectorCreate(MPI_COMM_WORLD, amg->global_size, amg->row_starts);
  amg->x   = hypre_ParVectorCreate(MPI_COMM_WORLD, amg->global_size, amg->row_starts);

  //point amg2013 to these assembled vectors
  xLocal   = hypre_ParVectorLocalVector(amg->x);
  rhsLocal = hypre_ParVectorLocalVector(amg->rhs);
  hypre_VectorData(xLocal)    = amg->xAssembled;
  hypre_VectorData(rhsLocal) = amg->rhsAssembled;

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
  
  amg->A = (HYPRE_ParCSRMatrix) A;
  amg->par_x = (HYPRE_ParVector) amg->x;
  amg->par_rhs = (HYPRE_ParVector) amg->rhs;


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
             (HYPRE_Vector) amg->par_rhs, (HYPRE_Vector) amg->par_x );

  return (void*) amg;
}


void * amg2013Setup(int    Nnum,            //unassembled size (Nverts*Nelements)
                    int    *row_starts,     //[numproc+1] global partition array 
                    int    nnz, 
                    int    *Ai,      // coo sparse matrix (globally indexed)
                    int    *Aj,
                    void   *Avals,
                    int    *sendSortId, 
                    int    *globalSortId, 
                    int    *compressId,
                    int    *sendCounts, 
                    int    *sendOffsets, 
                    int    *recvCounts, 
                    int    *recvOffsets,
                    const char* iintType, 
                    const char* dfloatType) {

  void *amg;

  int num_procs, myid;
  int *diag_i, *diag_j;
  double *diag_data;
  int *offd_i; 
  int *colMap=NULL;
  double *offd_data=NULL;
  int n;

  if (!strcmp(iintType,"long")) { 
    printf("Exiting due to use of ulong in iintType %s\n", iintType);
    exit(-1);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  int localStart = row_starts[myid];
  int localEnd   = row_starts[myid+1]-1;
  int numLocalRows = row_starts[myid+1]-row_starts[myid];

  diag_i = (int*) calloc(numLocalRows+1,sizeof(int));
  offd_i = (int*) calloc(numLocalRows+1,sizeof(int));

  //count local and nonlocal entries
  int cnt   = 0;
  int o_cnt = 0;
  int row   = 1;
  for (n=0;n<nnz;n++) { 
    if (Ai[n]>=row+localStart) { //next row 
      diag_i[row+1] = diag_i[row];
      offd_i[row+1] = offd_i[row];
      row++;
      n--;
    } else {
      if ((Aj[n]<localStart)||(Aj[n]>localEnd)){ //nonlocal column
        offd_i[row]++;
        o_cnt++;
      } else { //local colunm
        diag_i[row]++;
        cnt++;
      }
    }
  }
  diag_i[numLocalRows] = cnt;
  offd_i[numLocalRows] = o_cnt;

  diag_j = (int*) calloc(cnt,sizeof(int));
  diag_data = (double *) calloc(cnt,sizeof(double));
  if (o_cnt) {
    colMap = (int*) calloc(o_cnt,sizeof(int));
    offd_data = (double *) calloc(o_cnt,sizeof(double));
  } 

  cnt = 0;
  o_cnt = 0;
  for (n=0;n<nnz;n++) { 
    if ((Aj[n]<localStart)||(Aj[n]>localEnd)){ //nonlocal column
      colMap[o_cnt] = Aj[n];
      if (!strcmp(dfloatType,"float")) 
        offd_data[o_cnt++] = ((float*) Avals)[n];
      else 
        offd_data[o_cnt++] = ((double*) Avals)[n];
    } else { //local colunm
      diag_j[cnt] = Aj[n] - localStart;
      if (!strcmp(dfloatType,"float")) 
        diag_data[cnt++] = ((float*) Avals)[n];
      else 
        diag_data[cnt++] = ((double*) Avals)[n];
    }
  }

  amg = amg2013SetupCSR(Nnum, row_starts,    
                        diag_i, diag_j, diag_data,
                        offd_i, colMap, offd_data,
                        sendSortId, globalSortId, compressId,
                        sendCounts, sendOffsets, 
                        recvCounts, recvOffsets,
                        iintType, dfloatType);
  
  free(diag_i);
  free(offd_i);
  free(diag_j);
  free(colMap);
  free(diag_data);
  free(offd_data);

  return amg; 
}



int amg2013Solve(void* x,
                 void* AMG,
                 void* rhs) {
  
  int    num_iterations;
  double final_res_norm;

  int n,id;


  amg2013_t *amg = (amg2013_t *) AMG;


  if (!strcmp(amg->dfloatType,"float")) { 
    for (n=0;n<amg->Nnum;n++) 
      amg->rhsUnassembled[n] = ((float*) rhs)[n];
  } else { //double
    for (n=0;n<amg->Nnum;n++) 
      amg->rhsUnassembled[n] = ((double*) rhs)[n];
  }

  //sort by owner
  for (n=0;n<amg->Nnum;n++) 
    amg->rhsSort[n] = amg->rhsUnassembled[amg->sendSortId[n]];

  //Scatter nodes to their owners
  MPI_Alltoallv(amg->rhsSort, amg->sendCounts, amg->sendOffsets, MPI_CHAR,
                amg->rhsUnassembled, amg->recvCounts, amg->recvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort by globalid 
  for (n=0;n<amg->recvNnum;n++) 
    amg->rhsSort[n] = amg->rhsUnassembled[amg->globalSortId[n]];

  //gather
  for(n=0;n<amg->numLocalRows;++n){
    amg->rhsAssembled[n] = 0.;
    amg->xAssembled[n] = 0.;
    for (id=amg->compressId[n];id<amg->compressId[n+1];id++) 
      amg->rhsAssembled[n] += amg->rhsSort[id];
  }

  //solve
  HYPRE_PCGSolve(amg->solver, (HYPRE_Matrix) amg->A,
             (HYPRE_Vector) amg->par_rhs, (HYPRE_Vector) amg->par_x);
  
  //get solve info (could be useful)
  HYPRE_PCGGetNumIterations( amg->solver, &num_iterations );
  HYPRE_PCGGetFinalRelativeResidualNorm( amg->solver, &final_res_norm );

  printf("num it = %d, final res = %g\n", num_iterations, final_res_norm);
  
  //scatter
  for(n=0;n<amg->numLocalRows;++n){
    for (id = amg->compressId[n];id<amg->compressId[n+1];id++) 
      amg->xSort[id] = amg->xAssembled[n]; 
  }

  //sort by original rank
  for (n=0;n<amg->Nnum;n++) 
    amg->xUnassembled[amg->globalSortId[n]] = amg->xSort[n];

  //Scatter nodes back to their original rank
  MPI_Alltoallv(amg->xUnassembled, amg->recvCounts, amg->recvOffsets, MPI_CHAR,
                amg->xSort, amg->sendCounts, amg->sendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort back to original ordering
  for (n=0;n<amg->Nnum;n++) 
    amg->xUnassembled[amg->sendSortId[n]] = amg->xSort[n];

  if (!strcmp(amg->dfloatType,"float")) { 
    for (n=0;n<amg->Nnum;n++) 
      ((float*) x)[n] = amg->xUnassembled[n];
  } else { //double
    for (n=0;n<amg->Nnum;n++)
      ((double*) x)[n] = amg->xUnassembled[n];
  }

  return 0;
}


int amg2013Free(void* AMG) {
  amg2013_t *amg = (amg2013_t *) AMG;

  HYPRE_ParCSRPCGDestroy(amg->solver);

  HYPRE_BoomerAMGDestroy(amg->precond);

  free(amg->row_starts);

  free(amg->diag_i);
  free(amg->diag_j);
  free(amg->diag_data);

  free(amg->offd_i);
  if (amg->offd_nnz) {
    free(amg->offd_j);
    free(amg->colMap);
    free(amg->offd_data);
  }

  HYPRE_ParCSRMatrixDestroy(amg->A);
  HYPRE_ParVectorDestroy(amg->par_rhs);
  HYPRE_ParVectorDestroy(amg->par_x); 

  free(amg->xAssembled);
  free(amg->rhsAssembled);

  free(amg->xUnassembled);
  free(amg->rhsUnassembled);
  free(amg->xSort);
  free(amg->rhsSort);

  free(amg->sendCounts);
  free(amg->recvCounts);
  free(amg->sendOffsets);
  free(amg->recvOffsets);
  free(amg->sendSortId);
  free(amg->globalSortId);
  free(amg->compressId);

  free(amg);

  return 0;
}
