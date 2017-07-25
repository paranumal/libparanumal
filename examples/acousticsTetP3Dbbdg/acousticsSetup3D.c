#include "acoustics3D.h"

iint factorial(iint n) {
  iint retval = 1;
  for (iint i = n; i > 1; --i) retval *= i;
  return retval;
}

void acousticsSetup3D(mesh3D *mesh){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // set time step
  mesh->finalTime = 0.5;
  dfloat cfl = .4; // depends on the stability region size

  // errorStep
  mesh->errorStep = 10;

  // set penalty parameter
  mesh->Lambda2 = 0.5;

  dfloat *EtoDT = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
  dfloat hmin = 1e9;
  for(iint e=0;e<mesh->Nelements;++e){ 
    EtoDT[e] = 1e9;  

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)

      dfloat hest = .5/(sJ*invJ);

      // dt ~ cfl (h/(N+1)^2)/(Lambda^2*fastest wave speed)
      dfloat dtEst = cfl*hest/((mesh->N[e]+1.)*(mesh->N[e]+1.)*mesh->Lambda2);

      hmin = mymin(hmin, hest);
      EtoDT[e] = mymin(EtoDT[e], dtEst);
    }
  }

  //use dt on each element to setup MRAB
  int maxLevels = 10;
  meshMRABSetupP3D(mesh,EtoDT,maxLevels);


  mesh->Nfields = 4;
  
  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->NpMax*mesh->Nfields, sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(3*mesh->Nelements*mesh->NpMax*mesh->Nfields,
				sizeof(dfloat));
  mesh->fQM = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields,
            sizeof(dfloat));
  mesh->fQP = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields,
            sizeof(dfloat));

  // fix this later (initial conditions)
  dfloat time = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    iint N = mesh->N[e];
    for(iint n=0;n<mesh->Np[N];++n){
      dfloat x = mesh->x[n + mesh->NpMax*e];
      dfloat y = mesh->y[n + mesh->NpMax*e];
      dfloat z = mesh->z[n + mesh->NpMax*e];
      
      iint cnt = e*mesh->NpMax*mesh->Nfields + n*mesh->Nfields;
      //acousticsGaussianPulse3D(x, y, z, time,
			//	mesh->q+cnt,
			//	mesh->q+cnt+1,
			//	mesh->q+cnt+2,
			//	mesh->q+cnt+3);
      mesh->q[cnt+0] = 0.;
      mesh->q[cnt+1] = 0.;
      mesh->q[cnt+2] = 0.;
      mesh->q[cnt+3] = 0.;
    }
  }

  //Transform to BB modal space
  dfloat qtmp[mesh->Nfields*mesh->NpMax];
  for (iint e =0;e<mesh->Nelements;e++){
    iint cnt = e*mesh->NpMax*mesh->Nfields;
    iint N = mesh->N[e];

    for (iint n=0; n<mesh->Np[N]; n++){
      qtmp[n*mesh->Nfields + 0] = mesh->q[cnt+n*mesh->Nfields+0];
      qtmp[n*mesh->Nfields + 1] = mesh->q[cnt+n*mesh->Nfields+1];
      qtmp[n*mesh->Nfields + 2] = mesh->q[cnt+n*mesh->Nfields+2];
      qtmp[n*mesh->Nfields + 3] = mesh->q[cnt+n*mesh->Nfields+3];
      mesh->q[cnt+n*mesh->Nfields+0] = 0.0;
      mesh->q[cnt+n*mesh->Nfields+1] = 0.0;
      mesh->q[cnt+n*mesh->Nfields+2] = 0.0;
      mesh->q[cnt+n*mesh->Nfields+3] = 0.0;
    }
    for (iint n=0;n<mesh->Np[N];n++){
      for (iint m=0; m<mesh->Np[N]; m++){
        mesh->q[cnt+n*mesh->Nfields + 0] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+0];
        mesh->q[cnt+n*mesh->Nfields + 1] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+1];
        mesh->q[cnt+n*mesh->Nfields + 2] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+2];
        mesh->q[cnt+n*mesh->Nfields + 3] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+3];
      }
    }
  }

  printf("hmin = %g\n", hmin);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", mesh->dt);

  // output mesh
  meshVTU3D(mesh, "foo.vtu");

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  
  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);
  //  sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  // sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = Serial");

  occa::kernelInfo kernelInfo;

  mesh->device.setup(deviceConfig);

  mesh->o_NelList = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  
  mesh->o_D0ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_D1ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_D2ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_D3ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_Dvals = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));

  mesh->o_L0ids  = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_L0vals = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_ELids  = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_ELvals = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));

  mesh->o_BBLower     = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_BBRaiseids  = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_BBRaiseVals = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));

  mesh->o_cubInterpT = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_cubProjectT = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_cubDrWT = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_cubDsWT = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_cubDtWT = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));

  iint NMax = mesh->NMax;

  mesh->VBplot = (dfloat**) malloc((NMax+1)*sizeof(dfloat*));

  for (iint nn=1; nn <= NMax; nn++) {
    // deriv operators: transpose from row major to column major
    iint *D0ids = (iint*) calloc(mesh->Np[nn]*4,sizeof(iint));
    iint *D1ids = (iint*) calloc(mesh->Np[nn]*4,sizeof(iint));
    iint *D2ids = (iint*) calloc(mesh->Np[nn]*4,sizeof(iint));
    iint *D3ids = (iint*) calloc(mesh->Np[nn]*4,sizeof(iint));
    dfloat *Dvals = (dfloat*) calloc(mesh->Np[nn]*4,sizeof(dfloat));  

    iint    *L0ids = (iint*)   calloc(mesh->Nfp[nn]*7,sizeof(iint));
    dfloat *L0vals = (dfloat*) calloc(mesh->Nfp[nn]*7,sizeof(dfloat)); // tridiag
    iint    *ELids = (iint*)   calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn],sizeof(iint));
    dfloat *ELvals = (dfloat*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn],sizeof(dfloat));
    
    for (iint i = 0; i < mesh->Np[nn]; ++i){
      for (iint j = 0; j < 4; ++j){
        D0ids[i+j*mesh->Np[nn]] = mesh->D0ids[nn][j+i*4];
        D1ids[i+j*mesh->Np[nn]] = mesh->D1ids[nn][j+i*4];
        D2ids[i+j*mesh->Np[nn]] = mesh->D2ids[nn][j+i*4];
        D3ids[i+j*mesh->Np[nn]] = mesh->D3ids[nn][j+i*4];      
        Dvals[i+j*mesh->Np[nn]] = mesh->Dvals[nn][j+i*4];    
      }
    }

    for (iint i = 0; i < mesh->Nfp[nn]; ++i){
      for (iint j = 0; j < 7; ++j){
         L0ids [i+j*mesh->Nfp[nn]] = mesh->L0ids [nn][j+i*7];
         L0vals[i+j*mesh->Nfp[nn]] = mesh->L0vals[nn][j+i*7];
      }
    }
    
    for (iint i = 0; i < mesh->Np[nn]; ++i){
      for (iint j = 0; j < mesh->max_EL_nnz[nn]; ++j){
        ELids [i + j*mesh->Np[nn]] = mesh->ELids [nn][j+i*mesh->max_EL_nnz[nn]];
        ELvals[i + j*mesh->Np[nn]] = mesh->ELvals[nn][j+i*mesh->max_EL_nnz[nn]];
      }
    }
    
    //Build Vandermond matrix for conversion to nodal basis for plotting
    mesh->VBplot[nn] = (dfloat*) malloc(mesh->Np[nn]*mesh->NpMax*sizeof(dfloat));
    for (iint n=0;n<mesh->NpMax;n++) {
      dfloat r = mesh->r[NMax][n];
      dfloat s = mesh->s[NMax][n];
      dfloat t = mesh->t[NMax][n];

      dfloat l0 = -0.5*(1.+r+s+t); dfloat l1 = 0.5*(1.+r); dfloat l2 = 0.5*(1.+s); dfloat l3 = 0.5*(1.+t);
      
      iint cnt = 0;
      for (iint i=0;i<=nn;i++){
        for (iint j=0;j<=nn-i;j++){
          for (iint k=0;k<=nn-i-j;k++){
            mesh->VBplot[nn][n*mesh->Np[nn]+cnt] = ((dfloat) factorial(nn)/(factorial(i)*factorial(j)
                                            *factorial(k)*factorial(nn-i-j-k)))
                                            *pow(l0,nn-i-j-k)*pow(l1,k)*pow(l2,j)*pow(l3,i);
            cnt++;
          }
        }
      }
    }

    //Change cubature Interp and Project matrices
    for (iint n=0;n<mesh->cubNp[nn];n++) {
      dfloat r = mesh->cubr[nn][n];
      dfloat s = mesh->cubs[nn][n];
      dfloat t = mesh->cubt[nn][n];

      dfloat l0 = -0.5*(1.+r+s+t); dfloat l1 = 0.5*(1.+r); dfloat l2 = 0.5*(1.+s); dfloat l3 = 0.5*(1.+t);
      
      iint cnt = 0;
      for (iint i=0;i<=nn;i++){
        for (iint j=0;j<=nn-i;j++){
          for (iint k=0;k<=nn-i-j;k++){
            mesh->cubInterp[nn][n*mesh->Np[nn]+cnt] = ((dfloat) factorial(nn)/(factorial(i)*factorial(j)
                                            *factorial(k)*factorial(nn-i-j-k)))
                                            *pow(l0,nn-i-j-k)*pow(l1,k)*pow(l2,j)*pow(l3,i);
            cnt++;
          }
        }
      }
    }

    dfloat S[mesh->Np[nn]*mesh->cubNp[nn]];
    for (iint n=0;n<mesh->Np[nn];n++) {
      for (iint m =0;m<mesh->cubNp[nn];m++) {
        S[n*mesh->cubNp[nn] + m] = mesh->cubProject[nn][n*mesh->cubNp[nn] + m];
      }
    }
    for (iint n=0;n<mesh->Np[nn];n++) {
      for (iint m =0;m<mesh->cubNp[nn];m++) {
        mesh->cubProject[nn][n*mesh->cubNp[nn] + m] = 0.;
        for (iint i =0;i<mesh->Np[nn];i++)
          mesh->cubProject[nn][n*mesh->cubNp[nn]+m] 
            += mesh->invVB[nn][n*mesh->Np[nn] + i]*S[i*mesh->cubNp[nn]+m];
      }
    }

    // build volume cubature matrix transposes
    iint cubNpBlocked = mesh->Np[nn]*((mesh->cubNp[nn]+mesh->Np[nn]-1)/mesh->Np[nn]);
    dfloat *cubDrWT = (dfloat*) calloc(cubNpBlocked*mesh->Np[nn], sizeof(dfloat));
    dfloat *cubDsWT = (dfloat*) calloc(cubNpBlocked*mesh->Np[nn], sizeof(dfloat));
    dfloat *cubDtWT = (dfloat*) calloc(cubNpBlocked*mesh->Np[nn], sizeof(dfloat));
    dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    for(iint n=0;n<mesh->Np[nn];++n){
      for(iint m=0;m<mesh->cubNp[nn];++m){
        cubDrWT[n+m*mesh->Np[nn]] = mesh->cubDrW[nn][n*mesh->cubNp[nn]+m];
        cubDsWT[n+m*mesh->Np[nn]] = mesh->cubDsW[nn][n*mesh->cubNp[nn]+m];
        cubDtWT[n+m*mesh->Np[nn]] = mesh->cubDtW[nn][n*mesh->cubNp[nn]+m];

        cubProjectT[n+m*mesh->Np[nn]] = mesh->cubProject[nn][n*mesh->cubNp[nn]+m];
        cubInterpT[m+n*mesh->cubNp[nn]] = mesh->cubInterp[nn][m*mesh->Np[nn]+n];
        //      printf("%g @ ", cubInterpT[m+n*mesh->cubNp]);
      }
    }

    mesh->o_D0ids[nn] = mesh->device.malloc(mesh->Np[nn]*4*sizeof(iint),D0ids);
    mesh->o_D1ids[nn] = mesh->device.malloc(mesh->Np[nn]*4*sizeof(iint),D1ids);
    mesh->o_D2ids[nn] = mesh->device.malloc(mesh->Np[nn]*4*sizeof(iint),D2ids);
    mesh->o_D3ids[nn] = mesh->device.malloc(mesh->Np[nn]*4*sizeof(iint),D3ids);
    mesh->o_Dvals[nn] = mesh->device.malloc(mesh->Np[nn]*4*sizeof(dfloat),Dvals);

    mesh->o_L0ids [nn] = mesh->device.malloc(mesh->Nfp[nn]*7*sizeof(iint),L0ids);
    mesh->o_L0vals[nn] = mesh->device.malloc(mesh->Nfp[nn]*7*sizeof(dfloat),L0vals);
    mesh->o_ELids [nn] = mesh->device.malloc(mesh->Np[nn]*mesh->max_EL_nnz[nn]*sizeof(iint),ELids);
    mesh->o_ELvals[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->max_EL_nnz[nn]*sizeof(dfloat),ELvals);

    int Nfpp1 =  ((nn+2)*(nn+3))/2;
    int Nfpm1 =  ((nn)*(nn+1))/2;
    mesh->o_BBLower[nn]     = mesh->device.malloc(mesh->Nfp[nn]*Nfpm1*sizeof(dfloat),mesh->BBLower[nn]);
    mesh->o_BBRaiseids[nn]  = mesh->device.malloc(Nfpp1*3*sizeof(iint),mesh->BBRaiseids[nn]);
    mesh->o_BBRaiseVals[nn] = mesh->device.malloc(Nfpp1*3*sizeof(dfloat),mesh->BBRaiseVals[nn]);

    mesh->o_cubInterpT[nn]  = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubInterpT);
    mesh->o_cubProjectT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubProjectT);
    mesh->o_cubDrWT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubDrWT);
    mesh->o_cubDsWT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubDsWT);
    mesh->o_cubDtWT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubDtWT);

    
    free(D0ids); free(D1ids); free(D2ids); free(D3ids); free(Dvals);
    free(L0ids); free(L0vals); free(ELids); free(ELvals);
  }

  #if WADG
    // set heterogeneous c^2 for WADG
    mesh->c2 = (dfloat*) calloc(mesh->Nelements*mesh->cubNpMax,sizeof(dfloat));

    for(iint e=0;e<mesh->Nelements;++e){ /* for each element */
      
      iint id = e*mesh->Nverts+0;
      
      dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = mesh->EX[id+1];
      dfloat xe3 = mesh->EX[id+2];
      dfloat xe4 = mesh->EX[id+3];
      
      dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = mesh->EY[id+1];
      dfloat ye3 = mesh->EY[id+2];
      dfloat ye4 = mesh->EY[id+3];

      dfloat ze1 = mesh->EZ[id+0]; /* y-coordinates of vertices */
      dfloat ze2 = mesh->EZ[id+1];
      dfloat ze3 = mesh->EZ[id+2];
      dfloat ze4 = mesh->EZ[id+3];
      
      iint N = mesh->N[e];

      for(iint n=0;n<mesh->cubNp[N];++n){ /* for each node */
        
        // cubature node coordinates
        dfloat rn = mesh->cubr[N][n]; 
        dfloat sn = mesh->cubs[N][n];
        dfloat tn = mesh->cubt[N][n];

        /* physical coordinate of interpolation node */
        dfloat x = -0.5*(rn+sn+tn+1.)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4 ;
        dfloat y = -0.5*(rn+sn+tn+1.)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4 ;
        dfloat z = -0.5*(rn+sn+tn+1.)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4 ;
        
        // smoothly varying (sinusoidal) wavespeed
        //printf("M_PI = %f\n",M_PI);
        if (z<0.f) {
          mesh->c2[n + mesh->cubNpMax*e] = 0.2;//1.0 + 0.5*sin(M_PI*y);
        } else {
          mesh->c2[n + mesh->cubNpMax*e] = 1.0;
        }
      }
    }

    mesh->o_c2 = mesh->device.malloc(mesh->Nelements*mesh->cubNpMax*sizeof(dfloat),
         mesh->c2);
  #endif

  mesh->o_N = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint), mesh->N);  

  printf("Nverts = %d, Nfaces = %d\n",mesh->Nverts,mesh->Nfaces);
  if (mesh->Nverts==8){     // hardcoded for hexes

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements*mesh->NpMax*mesh->Nvgeo*sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->NfpMax*mesh->Nsgeo*sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements*mesh->NpMax*mesh->Nggeo*sizeof(dfloat),
        mesh->ggeo);

  }else if (mesh->Nverts==4){     // for tets

    mesh->o_vgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat),
                          mesh->vgeo);

    mesh->o_sgeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat),
                          mesh->sgeo);

    mesh->o_ggeo =
      mesh->device.malloc(mesh->Nelements*mesh->Nggeo*sizeof(dfloat),
        mesh->ggeo);

  }else{
    printf("Nverts = %d: unknown element type!\n",mesh->Nverts);
  }

  mesh->o_vmapM =
    mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint),
      mesh->vmapM);

  mesh->o_vmapP =
    mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint),
      mesh->vmapP);

  mesh->o_mapP  = 
    mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint), 
      mesh->mapP);

  mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint), mesh->EToE);
  mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint), mesh->EToF);

  mesh->o_EToB =
    mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),
      mesh->EToB);

  mesh->o_x =
    mesh->device.malloc(mesh->Nelements*mesh->NpMax*sizeof(dfloat), mesh->x);

  mesh->o_y =
    mesh->device.malloc(mesh->Nelements*mesh->NpMax*sizeof(dfloat), mesh->y);

  mesh->o_z =
    mesh->device.malloc(mesh->Nelements*mesh->NpMax*sizeof(dfloat), mesh->z);

  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*NpMax for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat));
  }

  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->NpMax*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(3*mesh->NpMax*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_fQM = 
    mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),mesh->fQM);
  mesh->o_fQP = 
    mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),mesh->fQP);

  //set up pml
  acousticsPmlSetup3D(mesh);

  //set up source injection
  acousticsSourceSetup3D(mesh,kernelInfo);

  mesh->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_MRABhaloIds    = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
  mesh->o_MRABelIdsP   = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory*));
  mesh->o_MRABhaloIdsP = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory*));

  for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
    if (mesh->MRABNelements[lev])
      mesh->o_MRABelementIds[lev] = mesh->device.malloc(mesh->MRABNelements[lev]*sizeof(iint),
         mesh->MRABelementIds[lev]);
    if (mesh->MRABNhaloElements[lev])
      mesh->o_MRABhaloIds[lev] = mesh->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(iint),
         mesh->MRABhaloIds[lev]);

    mesh->o_MRABelIdsP[lev]   = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
    mesh->o_MRABhaloIdsP[lev] = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
    for (int p=1;p<=mesh->NMax;p++) {
      if (mesh->MRABNelP[lev][p]) 
        mesh->o_MRABelIdsP[lev][p]   = mesh->device.malloc(mesh->MRABNelP[lev][p]*sizeof(iint),
         mesh->MRABelIdsP[lev][p]);
      if (mesh->MRABNhaloEleP[lev][p])
        mesh->o_MRABhaloIdsP[lev][p] = mesh->device.malloc(mesh->MRABNhaloEleP[lev][p]*sizeof(iint),
         mesh->MRABhaloIdsP[lev][p]);
    }
  }

  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_NMax", mesh->NMax);
  kernelInfo.addDefine("p_Nq", mesh->NMax+1);
  kernelInfo.addDefine("p_NpMax", mesh->NpMax);
  kernelInfo.addDefine("p_NfpMax", mesh->NfpMax);
  kernelInfo.addDefine("p_cubNpMax",mesh->cubNpMax);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfpMax", mesh->NfpMax*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  kernelInfo.addDefine("p_Nggeo", mesh->Nggeo);

  kernelInfo.addDefine("p_max_EL_nnzMax", mesh->max_EL_nnz[NMax]); // for Bernstein Bezier lift

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_NZID", NZID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);
  kernelInfo.addDefine("p_IHID", IHID);


  kernelInfo.addDefine("p_G00ID", G00ID);
  kernelInfo.addDefine("p_G01ID", G01ID);
  kernelInfo.addDefine("p_G02ID", G02ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_G12ID", G12ID);
  kernelInfo.addDefine("p_G22ID", G22ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);
  kernelInfo.addDefine("p_TXID", TXID);

  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);
  kernelInfo.addDefine("p_TYID", TYID);

  kernelInfo.addDefine("p_RZID", RZID);
  kernelInfo.addDefine("p_SZID", SZID);
  kernelInfo.addDefine("p_TZID", TZID);

  kernelInfo.addDefine("p_JWID", JWID);

  
  kernelInfo.addDefine("p_Lambda2", mesh->Lambda2);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
  }

  if(sizeof(iint)==4){
    kernelInfo.addDefine("iint","int");
  }
  if(sizeof(iint)==8){
    kernelInfo.addDefine("iint","long long int");
  }

  if(mesh->device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo.addCompilerFlag("--ftz=true");
    kernelInfo.addCompilerFlag("--prec-div=false");
    kernelInfo.addCompilerFlag("--prec-sqrt=false");
    kernelInfo.addCompilerFlag("--use_fast_math");
    kernelInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
  }

  mesh->volumeKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->surfaceKernel = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->updateKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->traceUpdateKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->pmlVolumeKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->pmlSurfaceKernel = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->pmlUpdateKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->pmlTraceUpdateKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));


  for (iint p=1;p<=mesh->NMax;p++) {
    occa::kernelInfo newInfo = kernelInfo;

    newInfo.addDefine("p_N", p);
    newInfo.addDefine("p_Np", mesh->Np[p]);
    newInfo.addDefine("p_Nfp", mesh->Nfp[p]);
    int Nfpp1 = (p+2)*(p+3)/2;
    int Nfpm1 = (p)*(p+1)/2;
    newInfo.addDefine("p_Nfpp1", Nfpp1);
    newInfo.addDefine("p_Nfpm1", Nfpm1);
    newInfo.addDefine("p_cubNp", mesh->cubNp[p]);

    newInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz[p]);

    int maxNodes = mymax(mesh->Np[p], mesh->Nfp[p]*mesh->Nfaces);
    newInfo.addDefine("p_maxNodes", maxNodes);

    int NblockV = 512/mesh->Np[p]; // works for CUDA
    newInfo.addDefine("p_NblockV", NblockV);

    int NblockS = 512/maxNodes; // works for CUDA
    newInfo.addDefine("p_NblockS", NblockS);

    int NblockCub = 512/mesh->cubNp[p]; // works for CUDA
    newInfo.addDefine("p_NblockCub", NblockCub);

    int maxCubNodes = mymax(mesh->cubNp[p], maxNodes);
    newInfo.addDefine("p_maxCubNodes", maxCubNodes);

    mesh->volumeKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABVolumeP3D.okl",
                 "acousticsbbdgMRABVolumeP3D",
                 newInfo);

    mesh->surfaceKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABSurfaceP3D.okl",
                 "acousticsbbdgMRABSurfaceP3D",
                 newInfo);

    #if WADG
      mesh->updateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP3D.okl",
                 "acousticsMRABUpdateP3D_wadg",
                 newInfo);
      mesh->traceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP3D.okl",
                 "acousticsMRABTraceUpdateP3D_wadg",
                 newInfo);
      mesh->pmlUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdateP3D.okl",
               "acousticsMRABPmlUpdateP3D_wadg",
                 newInfo);
      mesh->pmlTraceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdateP3D.okl",
               "acousticsMRABPmlTraceUpdateP3D_wadg",
                 newInfo);
    #else 
      mesh->updateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP3D.okl",
                 "acousticsMRABUpdateP3D",
                 newInfo);
      mesh->traceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP3D.okl",
                 "acousticsMRABTraceUpdateP3D",
                 newInfo);
      mesh->pmlUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdateP3D.okl",
               "acousticsMRABPmlUpdateP3D",
                 newInfo);
      mesh->pmlTraceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABPmlUpdateP3D.okl",
               "acousticsMRABPmlTraceUpdateP3D",
                 newInfo);
    #endif

    mesh->pmlVolumeKernel[p] =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABPmlVolume3D.okl",
               "acousticsbbdgMRABPmlVolume3D",
               newInfo);
    mesh->pmlSurfaceKernel[p] =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABPmlSurface3D.okl",
               "acousticsbbdgMRABPmlSurface3D",
               newInfo);
  }
  
  mesh->haloExtractKernel =
    mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);
}
