
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshLoadReferenceNodesTriP2D(mesh2D *mesh, int N){

  mesh->Np  = (iint*) malloc((N+1)*sizeof(iint));
  mesh->Nfp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->plotNp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->plotNelements = (iint*) malloc((N+1)*sizeof(iint));
  mesh->cubNp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->intNfp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->max_EL_nnz = (iint*) malloc((N+1)*sizeof(iint));

  mesh->r  = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->s  = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->Dr = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->Ds = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->MM = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));

  mesh->faceNodes = (iint**) malloc((N+1)*(sizeof(iint*)));  
  mesh->LIFT = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));

  mesh->plotR      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotS      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotInterp = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotEToV   = (iint**)   malloc((N+1)*sizeof(iint*));

  mesh->cubr      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubs      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubw      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubInterp = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubDrW    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubDsW    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubProject= (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->intInterp = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->intLIFT   = (dfloat**) malloc((N+1)*sizeof(dfloat*));

  mesh->VB     = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->invVB  = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->D1ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->D2ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->D3ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->Dvals  = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->VBq    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->PBq    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->L0vals = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->ELids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->ELvals = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->BBLower     = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->BBRaiseids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->BBRaiseVals = (dfloat**) malloc((N+1)*sizeof(dfloat*));

  for (int nn=1;nn<=N;nn++) {

    char fname[BUFSIZ];
    sprintf(fname, DHOLMES "/nodes/triangleN%02d.dat", nn);

    FILE *fp = fopen(fname, "r");

    printf("%s \n", fname);

    char buf[BUFSIZ];
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp);
    int Ncheck;
    sscanf(buf, "%d", &Ncheck);
    if(Ncheck != nn) printf("bu55er - wrong data file\n");
    //mesh->N[n] = Ncheck;
    mesh->Nfp[nn] = nn+1;

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp);
    int Npcheck;
    sscanf(buf, "%d", &Npcheck);
    mesh->Np[nn] = Npcheck;

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->r[nn] = (dfloat*) calloc(mesh->Np[nn], sizeof(dfloat));
    mesh->s[nn] = (dfloat*) calloc(mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat, mesh->r[nn]+n, mesh->s[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->Dr[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->Dr[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); // read comment

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->Ds[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->Ds[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); // read comment

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->MM[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->MM[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); // read comment

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->faceNodes[nn] = (iint*) calloc(mesh->Nfp[nn]*mesh->Nfaces, sizeof(iint));
    for(int n=0;n<mesh->Nfaces*mesh->Nfp[nn];++n){
      fscanf(fp, iintFormat, mesh->faceNodes[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->LIFT[nn] = (dfloat*) calloc(mesh->Nfp[nn]*mesh->Nfaces*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Nfaces*mesh->Nfp[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->LIFT[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);

    // read number of plot nodes
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); 
    sscanf(buf, iintFormat, mesh->plotNp+nn);

    // read plot node coordinates (hard code triangles)
    mesh->plotR[nn] = (dfloat*) calloc(mesh->plotNp[nn], sizeof(dfloat));
    mesh->plotS[nn] = (dfloat*) calloc(mesh->plotNp[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->plotNp[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat, mesh->plotR[nn]+n, mesh->plotS[nn]+n);
    }
    
    // read plot interpolation matrix
    mesh->plotInterp[nn] = (dfloat*) calloc(mesh->plotNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->plotNp[nn];++n){
      for(int m=0;m<mesh->Np[nn];++m){
        fscanf(fp, dfloatFormat, mesh->plotInterp[nn]+n*mesh->Np[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

    // read number of elements in plot node triangulation
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); 
    sscanf(buf, iintFormat, mesh->plotNelements+nn);

    // read number of vertices per plot element
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, iintFormat, &(mesh->plotNverts));
    
    // build and read in plot node triangulation
    mesh->plotEToV[nn] = (iint*) calloc(mesh->plotNelements[nn]*mesh->plotNverts, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->plotNelements[nn];++n){
      for(int m=0;m<mesh->plotNverts;++m){
        fscanf(fp, iintFormat, mesh->plotEToV[nn]+m + mesh->plotNverts*n);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

    // read number of volume cubature nodes
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); 
    sscanf(buf, iintFormat, mesh->cubNp+nn);

    // read cub nodes and weights
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->cubr[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    mesh->cubs[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    mesh->cubw[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    for(int n=0;n<mesh->cubNp[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat dfloatFormat, mesh->cubr[nn]+n, mesh->cubs[nn]+n, mesh->cubw[nn]+n);
    } 

    // read volume cubature interpolation matrix
    mesh->cubInterp[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment

    for(int n=0;n<mesh->cubNp[nn];++n){
      for(int m=0;m<mesh->Np[nn];++m){
        fscanf(fp, dfloatFormat, mesh->cubInterp[nn]+n*mesh->Np[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

    // read cubature weak 'r' differentiation matrix
    mesh->cubDrW[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn];++n){
      for(int m=0;m<mesh->cubNp[nn];++m){
        fscanf(fp, dfloatFormat, mesh->cubDrW[nn]+n*mesh->cubNp[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }
    // read cubature weak 's' differentiation matrix
    mesh->cubDsW[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn];++n){
      for(int m=0;m<mesh->cubNp[nn];++m){
        fscanf(fp, dfloatFormat, mesh->cubDsW[nn]+n*mesh->cubNp[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

      // read cubature projection matrix
    mesh->cubProject[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn];++n){
      for(int m=0;m<mesh->cubNp[nn];++m){
        fscanf(fp, dfloatFormat, mesh->cubProject[nn]+n*mesh->cubNp[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }


    // read number of surface integration nodes
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); 
    sscanf(buf, iintFormat, mesh->intNfp+nn);

    // read surface intergration node interpolation matrix
    mesh->intInterp[nn] 
      = (dfloat*) calloc(mesh->intNfp[nn]*mesh->Nfaces*mesh->Nfp[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment

    for(int n=0;n<mesh->intNfp[nn]*mesh->Nfaces;++n){
      for(int m=0;m<mesh->Nfp[nn];++m){
        fscanf(fp, dfloatFormat, mesh->intInterp[nn]+n*mesh->Nfp[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

    // read lift matrix from surface integration to volume nodes
    mesh->intLIFT[nn] = (dfloat*) calloc(mesh->intNfp[nn]*mesh->Nfaces*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn];++n){
      for(int m=0;m<mesh->intNfp[nn]*mesh->Nfaces;++m){
        fscanf(fp, dfloatFormat, mesh->intLIFT[nn]+n*mesh->intNfp[nn]*mesh->Nfaces+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }

    // ======================= BB data [added by JC] ===========================
  #if 1
    // BB Vandermonde matrices
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->VB[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->VB[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->invVB[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->invVB[nn]+n);
    }

    // sparse barycentric differentiation matrices
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->D1ids[nn]= (iint*) calloc(mesh->Np[nn]*3, sizeof(iint));
    for (int n=0;n<mesh->Np[nn]*3;++n){    
      fscanf(fp, iintFormat, mesh->D1ids[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->D2ids[nn] = (iint*) calloc(mesh->Np[nn]*3, sizeof(iint));
    for (int n=0;n<mesh->Np[nn]*3;++n){
      fscanf(fp, iintFormat, mesh->D2ids[nn]+n);
    }
    
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->D3ids[nn] = (iint*) calloc(mesh->Np[nn]*3, sizeof(iint));
    for (int n=0;n<mesh->Np[nn]*3;++n){
      fscanf(fp, iintFormat, mesh->D3ids[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->Dvals[nn] = (dfloat*) calloc(mesh->Np[nn]*3, sizeof(dfloat));
    for (int n=0;n<mesh->Np[nn]*3;++n){
      fscanf(fp, dfloatFormat, mesh->Dvals[nn]+n);
    }
    
    // BB cubature projection matrices
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->VBq[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->cubNp[nn], sizeof(dfloat));
    for (int n=0;n<mesh->Np[nn]*mesh->cubNp[nn];++n){
      fscanf(fp, dfloatFormat, mesh->VBq[nn]+n);
    }
    
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->PBq[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->cubNp[nn], sizeof(dfloat));
    for (int n=0;n<mesh->Np[nn]*mesh->cubNp[nn];++n){
      fscanf(fp, dfloatFormat, mesh->PBq[nn]+n);
    }

    // BB L0 matrix values
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->L0vals[nn] = (dfloat*) calloc(mesh->Nfp[nn]*3, sizeof(dfloat));
    for (int n=0;n<mesh->Nfp[nn]*3;++n){
      fscanf(fp, dfloatFormat, mesh->L0vals[nn]+n);
    }

    // BB lift reduction (EL) matrix (sparse format)
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->max_EL_nnz[nn] = mesh->Nfp[nn] + 2;
    mesh->ELids[nn] = (iint*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn], sizeof(iint));
    for (int n=0;n<mesh->Np[nn]*mesh->max_EL_nnz[nn];++n){
      fscanf(fp, iintFormat, mesh->ELids[nn]+n);
    }
    
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->ELvals[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn], sizeof(dfloat));
    for (int n=0;n<mesh->Np[nn]*mesh->max_EL_nnz[nn];++n){
      fscanf(fp, dfloatFormat, mesh->ELvals[nn]+n);
    }
    
    // BB degree raise matrix (sparse format)
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    int Nfpp1 = nn+2;
    mesh->BBRaiseids[nn] = (iint*) calloc(Nfpp1*2, sizeof(iint));
    for (int n=0;n<Nfpp1*2;++n){
      fscanf(fp, iintFormat, mesh->BBRaiseids[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->BBRaiseVals[nn] = (dfloat*) calloc(Nfpp1*2, sizeof(dfloat));
    for (int n=0;n<Nfpp1*2;++n){
      fscanf(fp, dfloatFormat, mesh->BBRaiseVals[nn]+n);
    }

    //BB degree lower matrix
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    int Nfpm1 = nn;
    mesh->BBLower[nn] = (dfloat*) calloc(Nfpm1*mesh->Nfp[nn], sizeof(dfloat));
    for (int n=0;n<Nfpm1*mesh->Nfp[nn];++n){
      fscanf(fp, dfloatFormat, mesh->BBLower[nn]+n);
    }

    // ============ end BB stuff ==================
  #endif
    
  #if 0
    printf("Dr: \n");
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        printf("%g ", mesh->Dr[n*mesh->Np+m]);
      }
      printf("\n");
    }

    printf("Ds: \n");
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        printf("%g ", mesh->Ds[n*mesh->Np+m]);
      }
      printf("\n");
    }
  #endif
    

    fclose(fp);
  }
}




