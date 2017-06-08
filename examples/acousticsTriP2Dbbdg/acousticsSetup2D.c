#include "acoustics2D.h"

// NBN: toggle use of 2nd stream
#define USE_2_STREAMS
// #undef USE_2_STREAMS

iint factorial(iint n) {
  iint retval = 1;
  for (iint i = n; i > 1; --i) retval *= i;
  return retval;
}

void acousticsSetup2D(mesh2D *mesh){

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
  meshMRABSetupP2D(mesh,EtoDT,maxLevels);


  mesh->Nfields = 4;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->NpMax*mesh->Nfields,
            sizeof(dfloat));
  mesh->rhsq = (dfloat*) calloc(3*mesh->Nelements*mesh->NpMax*mesh->Nfields,
            sizeof(dfloat));
  mesh->fQM = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields,
            sizeof(dfloat));
  mesh->fQP = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields,
            sizeof(dfloat));

  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    iint N  = mesh->N[e];
    for(iint n=0;n<mesh->Np[N];++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->NpMax*e];
      dfloat y = mesh->y[n + mesh->NpMax*e];

      cnt = e*mesh->NpMax*mesh->Nfields + n*mesh->Nfields;
      acousticsGaussianPulse2D(x, y, t,
             mesh->q+cnt, mesh->q+cnt+1, mesh->q+cnt+2);
    }
  }

  //Transform to BB modal space
  dfloat qtmp[mesh->Nfields*mesh->NpMax];
  for (iint e =0;e<mesh->Nelements;e++){
    cnt = e*mesh->NpMax*mesh->Nfields;
    iint N = mesh->N[e];
    for (iint n=0; n<mesh->Np[N]; n++){
      qtmp[n*mesh->Nfields + 0] = mesh->q[cnt+n*mesh->Nfields+0];
      qtmp[n*mesh->Nfields + 1] = mesh->q[cnt+n*mesh->Nfields+1];
      qtmp[n*mesh->Nfields + 2] = mesh->q[cnt+n*mesh->Nfields+2];
      mesh->q[cnt+n*mesh->Nfields+0] = 0.0;
      mesh->q[cnt+n*mesh->Nfields+1] = 0.0;
      mesh->q[cnt+n*mesh->Nfields+2] = 0.0;
    }
    for (iint n=0;n<mesh->Np[N];n++){
      for (iint m=0; m<mesh->Np[N]; m++){
        mesh->q[cnt+n*mesh->Nfields + 0] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+0];
        mesh->q[cnt+n*mesh->Nfields + 1] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+1];
        mesh->q[cnt+n*mesh->Nfields + 2] += mesh->invVB[N][n*mesh->Np[N]+m]*qtmp[m*mesh->Nfields+2];
      }
    }
  }

  printf("hmin = %g\n", hmin);
  printf("cfl = %g\n", cfl);
  printf("dt = %g\n", mesh->dt);
  printf("max wave speed = %g\n", sqrt(3.)*mesh->sqrtRT);


  // OCCA build stuff
  char deviceConfig[BUFSIZ];

  // use rank to choose DEVICE
  sprintf(deviceConfig, "mode = CUDA, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = OpenCL, deviceID = 0, platformID = 1");
  //sprintf(deviceConfig, "mode = OpenMP, deviceID = %d", 0);
  //sprintf(deviceConfig, "mode = Serial");

  occa::kernelInfo kernelInfo;

  mesh->device.setup(deviceConfig);

  mesh->o_NelList = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  
  mesh->o_D1ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_D2ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_D3ids = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
  mesh->o_Dvals = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));

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

  iint NMax = mesh->NMax;

  mesh->VBplot = (dfloat**) malloc((NMax+1)*sizeof(dfloat*));

  for (iint nn=1; nn <= NMax; nn++) {
    // deriv operators: transpose from row major to column major
    iint *D1ids = (iint*) calloc(mesh->Np[nn]*3,sizeof(iint));
    iint *D2ids = (iint*) calloc(mesh->Np[nn]*3,sizeof(iint));
    iint *D3ids = (iint*) calloc(mesh->Np[nn]*3,sizeof(iint));
    dfloat *Dvals = (dfloat*) calloc(mesh->Np[nn]*3,sizeof(dfloat));  

    dfloat *L0vals = (dfloat*) calloc(mesh->Nfp[nn]*3,sizeof(dfloat)); // tridiag
    iint    *ELids =   (iint*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn],sizeof(iint));
    dfloat *ELvals = (dfloat*) calloc(mesh->Np[nn]*mesh->max_EL_nnz[nn],sizeof(dfloat));
    
    for (iint i = 0; i < mesh->Np[nn]; ++i){
      for (iint j = 0; j < 3; ++j){
        D1ids[i+j*mesh->Np[nn]] = mesh->D1ids[nn][j+i*3];
        D2ids[i+j*mesh->Np[nn]] = mesh->D2ids[nn][j+i*3];
        D3ids[i+j*mesh->Np[nn]] = mesh->D3ids[nn][j+i*3];      
        Dvals[i+j*mesh->Np[nn]] = mesh->Dvals[nn][j+i*3];    
      }
    }
    
    for (iint i = 0; i < mesh->Nfp[nn]; ++i){
      for (iint j = 0; j < 3; ++j){
         L0vals[i+j*mesh->Nfp[nn]] = mesh->L0vals[nn][j+i*3];
      }
    }
    
    for (iint i = 0; i < mesh->Np[nn]; ++i){
      for (iint j = 0; j < mesh->max_EL_nnz[nn]; ++j){
        ELids[i + j*mesh->Np[nn]] = mesh->ELids[nn][j+i*mesh->max_EL_nnz[nn]];
        ELvals[i + j*mesh->Np[nn]] = mesh->ELvals[nn][j+i*mesh->max_EL_nnz[nn]];
      }
    }
    
    //Build Vandermond matrix for conversion to nodal basis for plotting
    mesh->VBplot[nn] = (dfloat*) malloc(mesh->Np[nn]*mesh->NpMax*sizeof(dfloat));
    for (iint n=0;n<mesh->NpMax;n++) {
      dfloat r = mesh->r[NMax][n];
      dfloat s = mesh->s[NMax][n];

      dfloat l1 = -0.5*(r+s); dfloat l2 = 0.5*(1.0+r); dfloat l3 = 0.5*(1.0+s);
      
      iint cnt = 0;
      for (iint i=0;i<=nn;i++){
        for (iint j=0;j<=nn-i;j++){
          mesh->VBplot[nn][n*mesh->Np[nn]+cnt] = ((dfloat) factorial(nn)/(factorial(i)*factorial(j)*factorial(nn-i-j)))
                                            *pow(l1,nn-i-j)*pow(l2,j)*pow(l3,i);
          cnt++;
        }
      }
    }
    
    //Change cubature Interp and Project matrices
    for (iint n=0;n<mesh->cubNp[nn];n++) {
      dfloat r = mesh->cubr[nn][n];
      dfloat s = mesh->cubs[nn][n];

      dfloat l1 = -0.5*(r+s); dfloat l2 = 0.5*(1.0+r); dfloat l3 = 0.5*(1.0+s);
      
      iint cnt = 0;
      for (iint i=0;i<=nn;i++){
        for (iint j=0;j<=nn-i;j++){
          mesh->cubInterp[nn][n*mesh->Np[nn]+cnt] = ((dfloat) factorial(nn)/(factorial(i)*factorial(j)*factorial(nn-i-j)))
                                            *pow(l1,nn-i-j)*pow(l2,j)*pow(l3,i);
          cnt++;
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
    dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNp[nn]*mesh->Np[nn], sizeof(dfloat));
    for(iint n=0;n<mesh->Np[nn];++n){
      for(iint m=0;m<mesh->cubNp[nn];++m){
        cubDrWT[n+m*mesh->Np[nn]] = mesh->cubDrW[nn][n*mesh->cubNp[nn]+m];
        cubDsWT[n+m*mesh->Np[nn]] = mesh->cubDsW[nn][n*mesh->cubNp[nn]+m];


        cubProjectT[n+m*mesh->Np[nn]] = mesh->cubProject[nn][n*mesh->cubNp[nn]+m];
        cubInterpT[m+n*mesh->cubNp[nn]] = mesh->cubInterp[nn][m*mesh->Np[nn]+n];
        //      printf("%g @ ", cubInterpT[m+n*mesh->cubNp]);
      }
    }

    mesh->o_D1ids[nn] = mesh->device.malloc(mesh->Np[nn]*3*sizeof(iint),D1ids);
    mesh->o_D2ids[nn] = mesh->device.malloc(mesh->Np[nn]*3*sizeof(iint),D2ids);
    mesh->o_D3ids[nn] = mesh->device.malloc(mesh->Np[nn]*3*sizeof(iint),D3ids);
    mesh->o_Dvals[nn] = mesh->device.malloc(mesh->Np[nn]*3*sizeof(dfloat),Dvals);

    mesh->o_L0vals[nn] = mesh->device.malloc(mesh->Nfp[nn]*3*sizeof(dfloat),L0vals);
    mesh->o_ELids[nn]  = mesh->device.malloc(mesh->Np[nn]*mesh->max_EL_nnz[nn]*sizeof(iint),ELids);
    mesh->o_ELvals[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->max_EL_nnz[nn]*sizeof(dfloat),ELvals);

    int Nfpp1 = nn+2; 
    int Nfpm1 = nn;
    mesh->o_BBLower[nn]     = mesh->device.malloc(mesh->Nfp[nn]*Nfpm1*sizeof(dfloat),mesh->BBLower[nn]);
    mesh->o_BBRaiseids[nn]  = mesh->device.malloc(Nfpp1*2*sizeof(iint),mesh->BBRaiseids[nn]);
    mesh->o_BBRaiseVals[nn] = mesh->device.malloc(Nfpp1*2*sizeof(dfloat),mesh->BBRaiseVals[nn]);
    
    mesh->o_cubInterpT[nn]  = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubInterpT);
    mesh->o_cubProjectT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubProjectT);
    mesh->o_cubDrWT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubDrWT);
    mesh->o_cubDsWT[nn] = mesh->device.malloc(mesh->Np[nn]*mesh->cubNp[nn]*sizeof(dfloat), cubDsWT);


    free(D1ids); free(D2ids); free(D3ids); free(Dvals);
    free(L0vals); free(ELids); free(ELvals);
  }

  #if WADG
    // set heterogeneous c^2 for WADG
    mesh->c2 = (dfloat*) calloc(mesh->Nelements*mesh->cubNpMax,sizeof(dfloat));

    for(iint e=0;e<mesh->Nelements;++e){ /* for each element */
      
      iint id = e*mesh->Nverts+0;
      
      dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
      dfloat xe2 = mesh->EX[id+1];
      dfloat xe3 = mesh->EX[id+2];
      
      dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
      dfloat ye2 = mesh->EY[id+1];
      dfloat ye3 = mesh->EY[id+2];
      
      iint N = mesh->N[e];

      for(iint n=0;n<mesh->cubNp[N];++n){ /* for each node */
        
        // cubature node coordinates
        dfloat rn = mesh->cubr[N][n]; 
        dfloat sn = mesh->cubs[N][n];

        /* physical coordinate of interpolation node */
        dfloat x = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
        dfloat y = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
        
        // smoothly varying (sinusoidal) wavespeed
        //printf("M_PI = %f\n",M_PI);
        if (y<0.f) {
          mesh->c2[n + mesh->cubNpMax*e] = 0.2;//1.0 + 0.5*sin(M_PI*y);
        } else {
          mesh->c2[n + mesh->cubNpMax*e] = 1.0;
        }
      }
    }

    mesh->o_c2 = mesh->device.malloc(mesh->Nelements*mesh->cubNpMax*sizeof(dfloat),
         mesh->c2);
  #endif

  mesh->o_N = mesh->device.malloc((mesh->totalHaloPairs+mesh->Nelements)*sizeof(iint), mesh->N);  

  mesh->o_vgeo = mesh->device.malloc(mesh->Nelements*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);  
  mesh->o_sgeo = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*mesh->Nsgeo*sizeof(dfloat), mesh->sgeo);

  mesh->o_vmapM = mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint), mesh->vmapM);
  mesh->o_vmapP = mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint), mesh->vmapP);
  mesh->o_mapP  = mesh->device.malloc(mesh->Nelements*mesh->NfpMax*mesh->Nfaces*sizeof(iint), mesh->mapP);

  mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint), mesh->EToE);
  mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint), mesh->EToF);
  mesh->o_EToB = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint), mesh->EToB);

  mesh->o_x = mesh->device.malloc(mesh->Nelements*mesh->NpMax*sizeof(dfloat), mesh->x);
  mesh->o_y = mesh->device.malloc(mesh->Nelements*mesh->NpMax*sizeof(dfloat), mesh->y);

  if(mesh->totalHaloPairs>0){
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      mesh->device.malloc(mesh->totalHaloPairs*sizeof(iint), mesh->haloElementList);
    
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat));
  }

  // find elements that have all neighbors on this process
  iint *internalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));
  iint *notInternalElementIds = (iint*) calloc(mesh->Nelements, sizeof(iint));

  iint Ninterior = 0, NnotInterior = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    iint flag = 0;
    
    for(iint f=0;f<mesh->Nfaces;++f) if(mesh->EToP[e*mesh->Nfaces+f]!=-1) flag = 1;

    if(!flag) { 
      internalElementIds[Ninterior++] = e;
    } else {
      notInternalElementIds[NnotInterior++] = e;
    }
  }

  printf("NinteriorElements = %d, NnotInternalElements = %d\n", Ninterior, NnotInterior);
  
  mesh->NinternalElements = Ninterior;
  mesh->NnotInternalElements = NnotInterior;
  mesh->o_internalElementIds    = mesh->device.malloc(Ninterior*sizeof(iint), internalElementIds);
  if(NnotInterior>0) mesh->o_notInternalElementIds = mesh->device.malloc(NnotInterior*sizeof(iint), notInternalElementIds);
  
  // OCCA allocate device memory (remember to go back for halo)
  mesh->o_q =
    mesh->device.malloc(mesh->NpMax*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);
  mesh->o_rhsq =
    mesh->device.malloc(3*mesh->NpMax*mesh->Nelements*mesh->Nfields*sizeof(dfloat), mesh->rhsq);
  mesh->o_fQM = 
    mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),mesh->fQM);
  mesh->o_fQP = 
    mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->NfpMax*mesh->Nfaces*mesh->Nfields*sizeof(dfloat),mesh->fQP);

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
    for (int p=0;p<=mesh->NMax;p++) {
      if (mesh->MRABNelP[lev][p]) 
        mesh->o_MRABelIdsP[lev][p]   = mesh->device.malloc(mesh->MRABNelP[lev][p]*sizeof(iint),
         mesh->MRABelIdsP[lev][p]);
      if (mesh->MRABNhaloEleP[lev][p])
        mesh->o_MRABhaloIdsP[lev][p] = mesh->device.malloc(mesh->MRABNhaloEleP[lev][p]*sizeof(iint),
         mesh->MRABhaloIdsP[lev][p]);
    }
  }

  //-------------------------------------
  // NBN: 2 streams for async MPI updates
  // {Vol, Surf, update}  run on q[0]
  // {halo-get, copy} run on q[1]
  //-------------------------------------
  mesh->stream0 = mesh->device.getStream();
#ifdef USE_2_STREAMS
  mesh->stream1 = mesh->device.createStream();  // NBN: second stream
#else
  mesh->stream1 = mesh->stream0;                // NBN: stream1 == stream0
#endif
  mesh->device.setStream(mesh->stream0);
  //-------------------------------------
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_NMax", mesh->NMax);
  kernelInfo.addDefine("p_Nq", mesh->NMax+1);
  kernelInfo.addDefine("p_NpMax", mesh->NpMax);
  kernelInfo.addDefine("p_NfpMax", mesh->NfpMax);
  kernelInfo.addDefine("p_Nfaces", mesh->Nfaces);
  kernelInfo.addDefine("p_NfacesNfpMax", mesh->NfpMax*mesh->Nfaces);
  kernelInfo.addDefine("p_Nvgeo", mesh->Nvgeo);
  kernelInfo.addDefine("p_Nsgeo", mesh->Nsgeo);
  kernelInfo.addDefine("p_Nggeo", mesh->Nggeo);

  kernelInfo.addDefine("p_NXID", NXID);
  kernelInfo.addDefine("p_NYID", NYID);
  kernelInfo.addDefine("p_SJID", SJID);
  kernelInfo.addDefine("p_IJID", IJID);
  kernelInfo.addDefine("p_WSJID", WSJID);

  kernelInfo.addDefine("p_intNfp", mesh->intNfp[mesh->NMax]);
  kernelInfo.addDefine("p_intNfpNfaces", mesh->intNfp[mesh->NMax]*mesh->Nfaces);

  if(sizeof(dfloat)==4){
    kernelInfo.addDefine("dfloat","float");
    kernelInfo.addDefine("dfloat4","float4");
    kernelInfo.addDefine("dfloat8","float8");
  }
  if(sizeof(dfloat)==8){
    kernelInfo.addDefine("dfloat","double");
    kernelInfo.addDefine("dfloat4","double4");
    kernelInfo.addDefine("dfloat8","double8");
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


  kernelInfo.addDefine("p_G00ID", G00ID);
  kernelInfo.addDefine("p_G01ID", G01ID);
  kernelInfo.addDefine("p_G11ID", G11ID);
  kernelInfo.addDefine("p_GWJID", GWJID);


  kernelInfo.addDefine("p_RXID", RXID);
  kernelInfo.addDefine("p_SXID", SXID);
  kernelInfo.addDefine("p_RYID", RYID);
  kernelInfo.addDefine("p_SYID", SYID);
  kernelInfo.addDefine("p_JID" , JID);

  kernelInfo.addDefine("p_NXID" , NXID);
  kernelInfo.addDefine("p_NYID" , NYID);
  kernelInfo.addDefine("p_SJID" , SJID);
  kernelInfo.addDefine("p_IJID" , IJID);

  kernelInfo.addDefine("p_NMax",mesh->NMax);
  kernelInfo.addDefine("p_NpMax",mesh->NpMax);
  kernelInfo.addDefine("p_NfpMax",mesh->NfpMax);
  kernelInfo.addDefine("p_cubNpMax",mesh->cubNpMax);

  kernelInfo.addDefine("p_Lambda2", 0.5f);

  mesh->volumeKernel       = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->surfaceKernel      = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->updateKernel       = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));
  mesh->traceUpdateKernel  = (occa::kernel *) malloc((mesh->NMax+1)*sizeof(occa::kernel));

  for (iint p=1;p<=mesh->NMax;p++) {
    occa::kernelInfo newInfo = kernelInfo;

    newInfo.addDefine("p_N", p);
    newInfo.addDefine("p_Np", mesh->Np[p]);
    newInfo.addDefine("p_Nfp", mesh->Nfp[p]);
    newInfo.addDefine("p_Nfpp1", p+2);
    newInfo.addDefine("p_Nfpm1", p);
    newInfo.addDefine("p_cubNp", mesh->cubNp[p]);

    newInfo.addDefine("p_max_EL_nnz", mesh->max_EL_nnz[p]);

    int maxNodes = mymax(mesh->Np[p], mesh->Nfp[p]*mesh->Nfaces);
    newInfo.addDefine("p_maxNodes", maxNodes);

    int NblockV = 512/mesh->Np[p]; // works for CUDA
    newInfo.addDefine("p_NblockV", NblockV);

    int NblockS = 512/maxNodes; // works for CUDA
    newInfo.addDefine("p_NblockS", NblockS);

    int maxCubNodes = mymax(mesh->cubNp[p], maxNodes);
    newInfo.addDefine("p_maxCubNodes", maxCubNodes);

    mesh->volumeKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABVolumeP2D.okl",
                 "acousticsbbdgMRABVolumeP2D",
                 newInfo);

    mesh->surfaceKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsbbdgMRABSurfaceP2D.okl",
                 "acousticsbbdgMRABSurfaceP2D",
                 newInfo);

    #if WADG
      mesh->updateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP2D.okl",
                 "acousticsMRABUpdateP2D_wadg",
                 newInfo);
      mesh->traceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP2D.okl",
                 "acousticsMRABTraceUpdateP2D_wadg",
                 newInfo);
    #else 
      mesh->updateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP2D.okl",
                 "acousticsMRABUpdateP2D",
                 newInfo);
      mesh->traceUpdateKernel[p] =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/acousticsMRABUpdateP2D.okl",
                 "acousticsMRABTraceUpdateP2D",
                 newInfo);
    #endif

  }
  
  mesh->haloExtractKernel =
      mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract2D.okl",
               "meshHaloExtract2D",
               kernelInfo);
}
