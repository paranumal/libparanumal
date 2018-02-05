#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/quad3dN%02d.dat", N);

  FILE *fp = fopen(fname, "r");
  
  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Ncheck;
  sscanf(buf, "%d", &Ncheck);
  if(Ncheck != N) printf("bu55er - wrong data file\n");
  mesh->N = N;
  mesh->Nq = N+1;
  mesh->Nfp = N+1;

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int Npcheck;
  sscanf(buf, "%d", &Npcheck);
  mesh->Np = Npcheck;

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->r = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  mesh->s = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->r+n, mesh->s+n);
  }

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-6;
  mesh->vertexNodes = (iint*) calloc(mesh->Nverts, sizeof(iint));
  for(iint n=0;n<mesh->Np;++n){
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[0] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]+1)*(mesh->s[n]+1)<NODETOL)
      mesh->vertexNodes[1] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[2] = n;
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]-1)*(mesh->s[n]-1)<NODETOL)
      mesh->vertexNodes[3] = n;
  }
  
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Dr = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Dr+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->Ds = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Np*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->Ds+n);
  }
  fgets(buf, BUFSIZ, fp); // read EOL

  fgets(buf, BUFSIZ, fp); // read comment
  mesh->faceNodes = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));
  for(int f=0;f<mesh->Nfaces;++f){
    for(int n=0;n<mesh->Nfp;++n){
      fscanf(fp, iintFormat, mesh->faceNodes+n + f*mesh->Nfp);
    }
  }
  fgets(buf, BUFSIZ, fp); // read EOL

    
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->LIFT = (dfloat*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->Nfaces*mesh->Nfp*mesh->Np;++n){
    fscanf(fp, dfloatFormat, mesh->LIFT+n);
  }
  fgets(buf, BUFSIZ, fp);

  /* 1D collocation differentiation matrix on GLL nodes */
  fgets(buf, BUFSIZ, fp); // read comment

  mesh->D = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    for(int m=0;m<mesh->N+1;++m){
      fscanf(fp, dfloatFormat, mesh->D+m+n*(mesh->N+1));
    }
  }
  fgets(buf, BUFSIZ, fp); // read comment

  /* 1D GLL node coordinates */
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->gllz = (dfloat*) calloc(mesh->N+1, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    fscanf(fp, dfloatFormat, mesh->gllz+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment
  
  /* 1D GLL node coordinates */
  fgets(buf, BUFSIZ, fp); // read comment
  mesh->gllw = (dfloat*) calloc(mesh->N+1, sizeof(dfloat));
  for(int n=0;n<mesh->N+1;++n){
    fscanf(fp, dfloatFormat, mesh->gllw+n);
  }
  fgets(buf, BUFSIZ, fp); // read comment

  // read number of plot nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->plotNp));

  // read plot node coordinates (hard code triangles)
  mesh->plotR = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  mesh->plotS = (dfloat*) calloc(mesh->plotNp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, dfloatFormat dfloatFormat, mesh->plotR+n, mesh->plotS+n);
  }
  
  // read plot interpolation matrix
  mesh->plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->plotInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read number of elements in plot node triangulation
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->plotNelements));

  // read number of vertices per plot element
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, iintFormat, &(mesh->plotNverts));
  
  // build and read in plot node triangulation
  mesh->plotEToV = (iint*) calloc(mesh->plotNelements*mesh->plotNverts, sizeof(iint));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->plotNelements;++n){
    for(int m=0;m<mesh->plotNverts;++m){
      fscanf(fp, iintFormat, mesh->plotEToV+m + mesh->plotNverts*n);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // projection info for OAS precon (one node overlap)
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, iintFormat, &(mesh->NqP));
  fgets(buf, BUFSIZ, fp);
  mesh->oasForward = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasForward+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf, BUFSIZ, fp);
  mesh->oasDiagOp = (dfloat*) calloc(mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    fscanf(fp, dfloatFormat, mesh->oasDiagOp+n);
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf,BUFSIZ,fp); // rest of line
  mesh->oasBack = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasBack+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // projection info for OAS precon (one node overlap)
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, iintFormat, &(mesh->NqP));
  fgets(buf, BUFSIZ, fp);
  mesh->oasForwardDg = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasForwardDg+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf, BUFSIZ, fp);
  mesh->oasDiagOpDg = (dfloat*) calloc(mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    fscanf(fp, dfloatFormat, mesh->oasDiagOpDg+n);
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  fgets(buf,BUFSIZ,fp); // rest of line
  mesh->oasBackDg = (dfloat*) calloc(mesh->NqP*mesh->NqP, sizeof(dfloat));
  for(int n=0;n<mesh->NqP;++n){
    for(int m=0;m<mesh->NqP;++m){
      fscanf(fp, dfloatFormat, mesh->oasBackDg+n*mesh->NqP+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  
  // read number of volume cubature nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->cubNp));

  // read volume cubature interpolation matrix
  mesh->cubInterp = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->cubNp;++n){
    for(int m=0;m<mesh->Np;++m){
      fscanf(fp, dfloatFormat, mesh->cubInterp+n*mesh->Np+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read cubature weak 'r' differentiation matrix
  mesh->cubDrW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDrW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }
  // read cubature weak 's' differentiation matrix
  mesh->cubDsW = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubDsW+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

    // read cubature projection matrix
  mesh->cubProject = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->cubNp;++m){
      fscanf(fp, dfloatFormat, mesh->cubProject+n*mesh->cubNp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }


  // read number of surface integration nodes
  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp); 
  sscanf(buf, iintFormat, &(mesh->intNfp));

  // read surface intergration node interpolation matrix
  mesh->intInterp 
    = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment
  for(int n=0;n<mesh->intNfp*mesh->Nfaces;++n){
    for(int m=0;m<mesh->Nfp;++m){
      fscanf(fp, dfloatFormat, mesh->intInterp+n*mesh->Nfp+m);
    }
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read lift matrix from surface integration to volume nodes
  mesh->intLIFT = (dfloat*) calloc(mesh->intNfp*mesh->Nfaces*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment

  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->intNfp*mesh->Nfaces;++m){
      fscanf(fp, dfloatFormat, mesh->intLIFT+n*mesh->intNfp*mesh->Nfaces+m);
      //      printf("%g ", mesh->intLIFT[n*mesh->intNfp*mesh->Nfaces+m]);
    }
    //    printf("\n");
    fgets(buf,BUFSIZ,fp); // rest of line
  }

  // read lift matrix from surface integration to volume nodes
  mesh->dualProjMatrix = (dfloat*) calloc(mesh->Nq*mesh->Nq*3, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp); // read comment

  for(int n=0;n<mesh->Nq;++n){
    for(int m=0;m<3*mesh->Nq;++m){
      fscanf(fp, dfloatFormat, mesh->dualProjMatrix+n*3*mesh->Nq+m);
      //      printf("%g ", mesh->dualProjMatrix[n*3*mesh->Nq+m]);
    }
    //    printf("\n");
    fgets(buf,BUFSIZ,fp); // rest of line
  }
  
  fclose(fp);
}
