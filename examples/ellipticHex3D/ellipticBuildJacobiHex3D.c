#include "ellipticHex3D.h"

void BuildLocalIpdgPatchAx(mesh3D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dlong eM, dfloat *A);

void BuildLocalContinuousPatchAx(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dfloat *A);

void ellipticBuildJacobiHex3D(solver_t* solver, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options){

  mesh3D *mesh = solver->mesh;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // build some monolithic basis arrays
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nk=0;nk<mesh->N+1;++nk){
    for(int nj=0;nj<mesh->N+1;++nj){
      for(int ni=0;ni<mesh->N+1;++ni){

        int node = 0;

        for(int k=0;k<mesh->N+1;++k){
          for(int j=0;j<mesh->N+1;++j){
            for(int i=0;i<mesh->N+1;++i){

              if(nk==k && nj==j && ni==i)
                B[mode*mesh->Np+node] = 1;
              if(nj==j && nk==k)
                Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
              if(ni==i && nk==k)
                Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
              if(ni==i && nj==j)
                Bt[mode*mesh->Np+node] = mesh->D[nk+mesh->Nq*k]; 
              
              ++node;
            }
          }
        }
        
        ++mode;
      }
    }
  }

  dlong diagNnum = mesh->Np*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  if(rank==0) printf("Building diagonal...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel
  {
    //temp patch storage
    dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

    #pragma omp for 
    for(dlong eM=0;eM<mesh->Nelements;++eM){
      //build the patch A matrix for this element
      if (strstr(options,"IPDG")) {
        BuildLocalIpdgPatchAx(mesh, tau, lambda, BCType, B, Br, Bs, Bt, eM, patchA);
      } else if (strstr(options,"CONTINUOUS")) {
        BuildLocalContinuousPatchAx(solver, lambda, eM, B, Br, Bs, Bt, patchA);
      }

      // compute the diagonal entries
      for(int n=0;n<mesh->Np;++n){
        diagA[eM*mesh->Np + n] = patchA[n*mesh->Np+n]; //store the inverse diagonal entry
      }
    }

    free(patchA);
  }

  if (strstr(options,"CONTINUOUS")) 
    gsParallelGatherScatter(mesh->hostGsh, diagA, dfloatString, "add"); 
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(rank==0) printf("done.\n");

  free(B); free(Br); free(Bs); free(Bt);
  free(diagA);
}
