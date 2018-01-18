
#include "ellipticTri2D.h"

void BuildLocalIpdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);
void BuildLocalBRdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);

void BuildLocalContinuousPatchAx(solver_t* solver, mesh2D* mesh, dfloat lambda,
                                  iint eM, dfloat *Srr, dfloat *Srs, 
                                  dfloat *Sss, dfloat *MM, dfloat *A);


void ellipticBuildJacobiTri2D(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options){

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Nfp;n++) {
      iint fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (iint m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (iint i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  dfloat *Srr, *Srs, *Sss, *MM;
  if (strstr(options,"CONTINUOUS")) {
    Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
    MM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

    if (strstr(options,"SPARSE")) {
      for (int k=0;k<mesh->SparseNnzPerRow;k++) {
        for (int n=0;n<mesh->Np;n++) {
          int m = mesh->sparseStackedNZ[n+k*mesh->Np]-1;
          if (m>-1) {
            Srr[m+n*mesh->Np] = mesh->sparseSrrT[n+k*mesh->Np];
            Srs[m+n*mesh->Np] = mesh->sparseSrsT[n+k*mesh->Np];
            Sss[m+n*mesh->Np] = mesh->sparseSssT[n+k*mesh->Np];  
          }
        }
      }
      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          MM[m+n*mesh->Np] = mesh->sparseMM[m+n*mesh->Np];
        }
      }
    } else {
      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          Srr[m+n*mesh->Np] = mesh->Srr[m+n*mesh->Np];
          Srs[m+n*mesh->Np] = mesh->Srs[m+n*mesh->Np] + mesh->Ssr[m+n*mesh->Np];
          Sss[m+n*mesh->Np] = mesh->Sss[m+n*mesh->Np];
          MM[m+n*mesh->Np] = mesh->MM[m+n*mesh->Np];
        }
      }
    }

  }

  iint diagNnum = basisNp*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){
    //build the patch A matrix for this element
    if (strstr(options,"IPDG")) {
      BuildLocalIpdgPatchAx(solver, mesh, basisNp, basis, tau, lambda, BCType, MS, eM, patchA);
    } else if (strstr(options,"BRDG")) {
      BuildLocalBRdgPatchAx(solver, mesh, basisNp, basis, tau, lambda, BCType, MS, eM, patchA);
    } else if (strstr(options,"CONTINUOUS")) {
      BuildLocalContinuousPatchAx(solver, mesh, lambda, eM, Srr, Srs, Sss, MM, patchA);
    }

    for(iint n=0;n<mesh->Np;++n) {
      diagA[eM*mesh->Np + n] = patchA[n*mesh->Np+n]; //store the diagonal entry
    }
  }

  if (strstr(options,"CONTINUOUS")) {
    if (strstr(options,"SPARSE")) for (iint n=0;n<mesh->Nelements*mesh->Np;n++) diagA[n] *= mesh->mapSgn[n];
    gsParallelGatherScatter(solver->hostGsh, diagA, dfloatString, "add"); 
    if (strstr(options,"SPARSE")) for (iint n=0;n<mesh->Nelements*mesh->Np;n++) diagA[n] *= mesh->mapSgn[n];
  }

  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (iint n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (strstr(options,"CONTINUOUS")&&(mesh->mask[n]==0)) continue;
    (*invDiagA)[n] = 1/diagA[n];
  }

  if (strstr(options,"CONTINUOUS")) {
    free(Srr);
    free(Srs);
    free(Sss);
    free(MM);
  }

  free(diagA);
  free(patchA);
  free(MS);
}
