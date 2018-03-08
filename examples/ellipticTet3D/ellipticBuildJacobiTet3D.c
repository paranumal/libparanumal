#include "ellipticTet3D.h"

void BuildLocalIpdgPatchAx(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A);

void BuildLocalContinuousPatchAx(solver_t* solver, mesh3D* mesh, dfloat lambda, dlong eM, 
                                  dfloat *Srr, dfloat *Srs, dfloat *Srt, 
                                  dfloat *Sss, dfloat *Sst, dfloat *Stt, 
                                  dfloat *MM, dfloat *A);

void ellipticBuildJacobiTet3D(solver_t *solver, mesh3D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (int m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  dfloat *Srr, *Srs, *Srt, *Sss, *Sst, *Stt, *MM;
  if (strstr(options,"CONTINUOUS")) {
    Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Srt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Sst = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Stt = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    MM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
        Srr[m+n*mesh->Np] = mesh->Srr[m+n*mesh->Np];
        Srs[m+n*mesh->Np] = mesh->Srs[m+n*mesh->Np] + mesh->Ssr[m+n*mesh->Np];
        Srt[m+n*mesh->Np] = mesh->Srt[m+n*mesh->Np] + mesh->Str[m+n*mesh->Np];
        Sss[m+n*mesh->Np] = mesh->Sss[m+n*mesh->Np];
        Sst[m+n*mesh->Np] = mesh->Sst[m+n*mesh->Np] + mesh->Sts[m+n*mesh->Np];
        Stt[m+n*mesh->Np] = mesh->Stt[m+n*mesh->Np];
        MM[m+n*mesh->Np] = mesh->MM[m+n*mesh->Np];
      }
    }
  }

  dlong diagNnum = basisNp*mesh->Nelements;

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
        BuildLocalIpdgPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, eM, patchA);
      } else if (strstr(options,"CONTINUOUS")) {
        BuildLocalContinuousPatchAx(solver, mesh, lambda, eM, Srr, Srs, Srt, Sss, Sst, Stt, MM, patchA);
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

  if (strstr(options,"CONTINUOUS")) {
    free(Srr);
    free(Srs);
    free(Srt);
    free(Sss);
    free(Sst);
    free(Stt);
    free(MM);
  }

  free(diagA);
  free(MS);
}
