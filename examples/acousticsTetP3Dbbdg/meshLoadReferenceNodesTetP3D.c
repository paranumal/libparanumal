
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshLoadReferenceNodesTetP3D(mesh3D *mesh, int N){

  //local element order
  mesh->N = (iint*) malloc(mesh->Nelements*sizeof(iint));

  for (int e=0;e<mesh->Nelements;e++){
    mesh->N[e] = N;
    if (e%2==0) mesh->N[e] = N-1;
  }

  mesh->Np  = (iint*) malloc((N+1)*sizeof(iint));
  mesh->Nfp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->plotNp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->plotNelements = (iint*) malloc((N+1)*sizeof(iint));
  mesh->cubNp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->intNfp = (iint*) malloc((N+1)*sizeof(iint));
  mesh->max_EL_nnz = (iint*) malloc((N+1)*sizeof(iint));

  mesh->r  = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->s  = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->t  = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->Dr = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->Ds = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));
  mesh->Dt = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));

  mesh->faceNodes = (iint**) malloc((N+1)*(sizeof(iint*)));  
  mesh->LIFT = (dfloat**) malloc((N+1)*(sizeof(dfloat*)));

  mesh->plotR      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotS      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotT      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotInterp = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->plotEToV   = (iint**)   malloc((N+1)*sizeof(iint*));

  mesh->cubr      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubs      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubt      = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubInterp = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubDrW    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubDsW    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubDtW    = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->cubProject= (dfloat**) malloc((N+1)*sizeof(dfloat*));

  mesh->VB     = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->invVB  = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->D0ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->D1ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->D2ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->D3ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->Dvals  = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->L0ids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->L0vals = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->ELids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->ELvals = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->BBLower     = (dfloat**) malloc((N+1)*sizeof(dfloat*));
  mesh->BBRaiseids  = (iint**)   malloc((N+1)*sizeof(iint*));
  mesh->BBRaiseVals = (dfloat**) malloc((N+1)*sizeof(dfloat*));

  for (int nn=1;nn<=N;nn++) {

    char fname[BUFSIZ];
    sprintf(fname, DHOLMES "/nodes/tetN%02d.dat", nn);

    FILE *fp = fopen(fname, "r");

    char buf[BUFSIZ];
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp);
    int Ncheck;
    sscanf(buf, "%d", &Ncheck);
    if(Ncheck != nn) printf("bu55er - wrong data file\n");

    //mesh->N = N;
    mesh->Np[nn] = ((nn+1)*(nn+2)*(nn+3))/6;
    mesh->Nfp[nn] = ((nn+1)*(nn+2))/2;

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp);
    int Npcheck;
    sscanf(buf, "%d", &Npcheck);
    //mesh->Np = Npcheck;

    fgets(buf, BUFSIZ, fp); // read comment
    mesh->r[nn] = (dfloat*) calloc(mesh->Np[nn], sizeof(dfloat));
    mesh->s[nn] = (dfloat*) calloc(mesh->Np[nn], sizeof(dfloat));
    mesh->t[nn] = (dfloat*) calloc(mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat dfloatFormat,
  	   mesh->r[nn]+n, mesh->s[nn]+n, mesh->t[nn]+n);
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
    mesh->Dt[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->Dt[nn]+n);
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
    fgets(buf, BUFSIZ, fp); // read comment
    
    // read number of plot nodes
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); 
    sscanf(buf, iintFormat, mesh->plotNp+nn);

    // read plot node coordinates (hard code triangles)
    mesh->plotR[nn] = (dfloat*) calloc(mesh->plotNp[nn], sizeof(dfloat));
    mesh->plotS[nn] = (dfloat*) calloc(mesh->plotNp[nn], sizeof(dfloat));
    mesh->plotT[nn] = (dfloat*) calloc(mesh->plotNp[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->plotNp[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat dfloatFormat, mesh->plotR[nn]+n, mesh->plotS[nn]+n, mesh->plotT[nn]+n);
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

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read number of cubature points
    sscanf(buf, iintFormat, mesh->cubNp+nn);
    printf("cubNp[%d] = %d \n", nn, mesh->cubNp[nn]);
      
    mesh->cubr[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    mesh->cubs[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    mesh->cubt[nn] = (dfloat*) calloc(mesh->cubNp[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->cubNp[nn];++n){
      fgets(buf, BUFSIZ, fp);
      sscanf(buf, dfloatFormat dfloatFormat dfloatFormat, mesh->cubr[nn]+n, mesh->cubs[nn]+n, mesh->cubt[nn]+n);
    }
    
    mesh->cubInterp[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->cubNp[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->cubInterp[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); // read comment  

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

     // read cubature weak 's' differentiation matrix
    mesh->cubDtW[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn];++n){
      for(int m=0;m<mesh->cubNp[nn];++m){
        fscanf(fp, dfloatFormat, mesh->cubDtW[nn]+n*mesh->cubNp[nn]+m);
      }
      fgets(buf,BUFSIZ,fp); // rest of line
    }
    
    mesh->cubProject[nn] = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->cubNp[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->cubProject[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); 

    //-------------Berstein Bezier DG stuff added by NC--------------------//
    mesh->VB[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->VB[nn]+n);
    }
    fgets(buf, BUFSIZ, fp); 
    
    mesh->invVB[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->Np[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*mesh->Np[nn];++n){
      fscanf(fp, dfloatFormat, mesh->invVB[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->D0ids[nn] = (iint*) calloc(mesh->Np[nn]*4, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*4;++n){
      fscanf(fp, iintFormat, mesh->D0ids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->D1ids[nn] = (iint*) calloc(mesh->Np[nn]*4, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*4;++n){
      fscanf(fp, iintFormat, mesh->D1ids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->D2ids[nn] = (iint*) calloc(mesh->Np[nn]*4, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*4;++n){
      fscanf(fp, iintFormat, mesh->D2ids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->D3ids[nn] = (iint*) calloc(mesh->Np[nn]*4, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*4;++n){
      fscanf(fp, iintFormat, mesh->D3ids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->Dvals[nn] = (dfloat*) calloc(mesh->Np[nn]*4, sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*4;++n){
      fscanf(fp, dfloatFormat, mesh->Dvals[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->L0ids[nn] = (iint*) calloc(mesh->Nfp[nn]*7, sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Nfp[nn]*7;++n){
      fscanf(fp, iintFormat, mesh->L0ids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->L0vals[nn] = (dfloat*) calloc(mesh->Nfp[nn]*7, sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Nfp[nn]*7;++n){
      fscanf(fp, dfloatFormat, mesh->L0vals[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->max_EL_nnz[nn] = mesh->Nfp[nn]+3;
    mesh->ELids[nn] = (iint*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn], sizeof(iint));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*mesh->max_EL_nnz[nn];++n){
      fscanf(fp, iintFormat, mesh->ELids[nn]+n);
    }
    fgets(buf, BUFSIZ, fp);   

    mesh->ELvals[nn] = (dfloat*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn], sizeof(dfloat));
    fgets(buf, BUFSIZ, fp); // read comment
    for(int n=0;n<mesh->Np[nn]*mesh->max_EL_nnz[nn];++n){
      fscanf(fp, dfloatFormat, mesh->ELvals[nn]+n);
    }

    // BB degree raise matrix (sparse format)
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->BBRaiseids[nn] = (iint*) calloc(mesh->Nfp[nn]*3, sizeof(iint));
    for (int n=0;n<mesh->Nfp[nn]*3;++n){
      fscanf(fp, iintFormat, mesh->BBRaiseids[nn]+n);
    }

    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    mesh->BBRaiseVals[nn] = (dfloat*) calloc(mesh->Nfp[nn]*3, sizeof(dfloat));
    for (int n=0;n<mesh->Nfp[nn]*3;++n){
      fscanf(fp, dfloatFormat, mesh->BBRaiseVals[nn]+n);
    }

    //BB degree lower matrix
    fgets(buf, BUFSIZ, fp); // read comment
    fgets(buf, BUFSIZ, fp); // read comment
    iint NfpPlusOne =  ((nn+2)*(nn+3))/2;
    mesh->BBLower[nn] = (dfloat*) calloc(NfpPlusOne*mesh->Nfp[nn], sizeof(dfloat));
    for (int n=0;n<NfpPlusOne*mesh->Nfp[nn];++n){
      fscanf(fp, dfloatFormat, mesh->BBLower[nn]+n);
    }
    
  #if 0
    for(int n=0;n<mesh->cubNp;++n){
      printf("rq,sq,tq = %f, %f, %f\n",mesh->cubr[n],mesh->cubs[n],mesh->cubt[n]);
    }

    printf("Vq: \n");
    for(int n=0;n<mesh->cubNp;++n){
      for(int m=0;m<mesh->Np;++m){
        printf("%f ",mesh->cubInterp[n*mesh->Np+m]);
      }
      printf("\n");
    }

    printf("Pq: \n");
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->cubNp;++m){
        printf("%f ",mesh->cubProject[n*mesh->cubNp+m]);
      }
      printf("\n");
    }
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

    printf("Dt: \n");
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        printf("%g ", mesh->Dt[n*mesh->Np+m]);
      }
      printf("\n");
    }
    
    printf("LIFT: \n");
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Nfaces*mesh->Nfp;++m){
        printf("%g ", mesh->LIFT[n*mesh->Nfaces*mesh->Nfp+m]);
      }
      printf("\n");
    }

    printf("faceNodes: \n");
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nfp;++n){
        printf("%d ", mesh->faceNodes[n+f*mesh->Nfp]);
      }
      printf("\n");
    }
  #endif

    fclose(fp);
  }

  mesh->NMax = N;
  mesh->NpMax = mesh->Np[N];
  mesh->NfpMax = mesh->Nfp[N];  
  mesh->cubNpMax = mesh->cubNp[N];  
}
