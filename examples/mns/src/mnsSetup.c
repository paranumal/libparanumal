#include "mns.h"
#include "omp.h"
#include <unistd.h>

mns_t *mnsSetup(mesh_t *mesh, setupAide options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  if (size==1) options.getArgs("DEVICE NUMBER" ,deviceID);

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",deviceID);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", deviceID, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }
  
  //set number of omp threads to use
  int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int Nthreads = Ncores/totalDevices;
  omp_set_num_threads(Nthreads);
  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    printf("Rank %d: Ncores = %d, Nthreads = %d\n", rank, Ncores, Nthreads);



  mns_t *mns = (mns_t*) calloc(1, sizeof(mns_t));
  mns->mesh = mesh;
  mns->options = options;

  options.getArgs("MESH DIMENSION", mns->dim);
  options.getArgs("ELEMENT TYPE", mns->elementType);

  mns->NVfields = (mns->dim==3) ? 3:2; //  Total Number of Velocity Fields
  mns->NTfields = (mns->dim==3) ? 4:3; // Total Velocity + Pressure

  mesh->Nfields = 1; 

  mns->g0 =  1.0;

  if (options.compareArgs("TIME INTEGRATOR", "ARK1")) {
    mns->Nstages = 1;
    int Nrk = 2;
    dfloat rkC[2] = {0.0, 1.0};

    dfloat erkA[2*2] ={0.0, 0.0,\
                       1.0, 0.0};
    dfloat irkA[2*2] ={0.0, 0.0,\
                       0.0, 1.0};
    
    dfloat prkA[2*2] ={0.0, 0.0,\
                       0.0, 1.0};
    dfloat prkB[2*2] ={0.0, 0.0,\
                       1.0, 0.0};                       

    mns->Nrk = Nrk;
    mns->rkC = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->erkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->irkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkB = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));

    memcpy(mns->rkC, rkC, mns->Nrk*sizeof(dfloat));
    memcpy(mns->erkA, erkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->irkA, irkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkA, prkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkB, prkB, mns->Nrk*mns->Nrk*sizeof(dfloat));

    mns->g0 =  1.0;
    mns->embeddedRKFlag = 0; //no embedded method
  } else if (options.compareArgs("TIME INTEGRATOR", "ARK2")) {
    mns->Nstages = 2;
    int Nrk = 3;

    dfloat gamma = (2.0-sqrt(2.0))/2.0;
    dfloat delta = 1.0 - 1.0/(2.0*gamma);

    dfloat rkC[3]    ={  0.0,   gamma,   1.0};

    dfloat erkA[3*3] ={  0.0,     0.0,   0.0,\
                       gamma,     0.0,   0.0,\
                       delta, 1-delta,   0.0};
    dfloat irkA[3*3] ={  0.0,     0.0,   0.0,\
                         0.0,   gamma,   0.0,\
                         0.0, 1-gamma, gamma};
    
    dfloat prkA[3*3] ={  0.0,     0.0,   0.0,\
                         0.0,   gamma,   0.0,\
                         0.5,     0.0,   0.5};

    dfloat prkB[3*3] = {0.0,     0.0,   0.0,\
                        1.0,     0.0,   0.0,\
                        1.0,     0.0,   0.0};

    mns->Nrk = Nrk;
    mns->rkC = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->erkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->irkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkB = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    

    memcpy(mns->rkC, rkC, mns->Nrk*sizeof(dfloat));
    memcpy(mns->erkA, erkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->irkA, irkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkA, prkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkB, prkB, mns->Nrk*mns->Nrk*sizeof(dfloat));
    
    mns->g0 =  1.0/gamma;
    mns->embeddedRKFlag = 0; //no embedded method
  } else if (options.compareArgs("TIME INTEGRATOR", "ARK3")) {
    mns->Nstages = 4;
    int Nrk = 4;

    dfloat erkA[4*4] ={                              0.0,                              0.0,                               0.0, 0.0,\
                         1767732205903.0/2027836641118.0,                              0.0,                               0.0, 0.0,\
                        5535828885825.0/10492691773637.0,  788022342437.0/10882634858940.0,                               0.0, 0.0,\
                        6485989280629.0/16251701735622.0, -4246266847089.0/9704473918619.0, 10755448449292.0/10357097424841.0, 0.0};
    dfloat erkB[4] = {1471266399579.0/7840856788654.0, \
                      -4482444167858.0/7529755066697.0, \
                      11266239266428.0/11593286722821.0, \
                      1767732205903.0/4055673282236.0};
    dfloat erkE[4] = {1471266399579.0/7840856788654.0 - 2756255671327.0/12835298489170.0,\
                      -4482444167858.0/7529755066697.0 - -10771552573575.0/22201958757719.0,\
                      11266239266428.0/11593286722821.0 - 9247589265047.0/10645013368117.0,\
                      1767732205903.0/4055673282236.0 - 2193209047091.0/5459859503100.0};

    dfloat irkA[4*4] ={                              0.0,                              0.0,                               0.0,                             0.0,\
                         1767732205903.0/4055673282236.0,  1767732205903.0/4055673282236.0,                               0.0,                             0.0,\
                        2746238789719.0/10658868560708.0,  -640167445237.0/6845629431997.0,   1767732205903.0/4055673282236.0,                             0.0,\
                         1471266399579.0/7840856788654.0, -4482444167858.0/7529755066697.0, 11266239266428.0/11593286722821.0, 1767732205903.0/4055673282236.0};
    dfloat irkB[4] = {1471266399579.0/7840856788654.0,\
                      -4482444167858.0/7529755066697.0,\
                      11266239266428.0/11593286722821.0,\
                      1767732205903.0/4055673282236.0};
    dfloat irkE[4] = {1471266399579.0/7840856788654.0 - 2756255671327.0/12835298489170.0,\
                      -4482444167858.0/7529755066697.0 - -10771552573575.0/22201958757719.0,
                      11266239266428.0/11593286722821.0 - 9247589265047.0/10645013368117.0,\
                      1767732205903.0/4055673282236.0 - 2193209047091.0/5459859503100.0};

    dfloat rkC[4] = {0.0, \
                    1767732205903.0/2027836641118.0, \
                    3.0/5.0, \
                    1.0};

    dfloat prkA[4*4] ={  0.0,                             0.0,       0.0,   0.0,\
                         0.0, 1767732205903.0/2027836641118.0,       0.0,   0.0,\
                         0.5,                             0.0,       0.5,   0.5};

    dfloat prkB[4*4] ={  0.0,                             0.0,       0.0,   0.0,\
                         0.0, 1767732205903.0/2027836641118.0,       0.0,   0.0,\
                         0.5,                             0.0,       0.5,   0.5};

    mns->Nrk = Nrk;
    mns->erkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->irkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkA = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->prkB = (dfloat*) calloc(mns->Nrk*mns->Nrk, sizeof(dfloat));
    mns->erkB = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->erkE = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->irkB = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->irkE = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));
    mns->rkC  = (dfloat*) calloc(mns->Nrk, sizeof(dfloat));


    memcpy(mns->erkA, erkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->irkA, irkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkA, prkA, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->prkB, prkB, mns->Nrk*mns->Nrk*sizeof(dfloat));
    memcpy(mns->erkB, erkB, mns->Nrk*sizeof(dfloat));
    memcpy(mns->erkE, erkE, mns->Nrk*sizeof(dfloat));
    memcpy(mns->irkB, irkB, mns->Nrk*sizeof(dfloat));
    memcpy(mns->irkE, irkE, mns->Nrk*sizeof(dfloat));
    memcpy(mns->rkC, rkC, mns->Nrk*sizeof(dfloat));
    
    mns->g0 =  4055673282236.0/1767732205903.0;
    mns->embeddedRKFlag = 1; 
  } 


  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    mns->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    mns->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    mns->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));

    mns->extC = (dfloat*) calloc(3, sizeof(dfloat));
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    mns->Nstages = 1;
    mns->temporalOrder = 1;
    mns->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    mns->Nstages = 2;
    mns->temporalOrder = 2;
    mns->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    mns->Nstages = 3;
    mns->temporalOrder = 3;
    mns->g0 = 11.f/6.f;
  }

  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  mns->Ntotal = Ntotal;
  mns->fieldOffset = Ntotal;
  mns->Nblock = (Nlocal+blockSize-1)/blockSize;

  // compute samples of q at interpolation nodes
  mns->U      = (dfloat*) calloc(mns->NVfields*mns->Nstages*Ntotal,sizeof(dfloat));
  mns->P      = (dfloat*) calloc(              mns->Nstages*Ntotal,sizeof(dfloat));
  mns->Phi    = (dfloat*) calloc(                           Ntotal,sizeof(dfloat));

  //rhs storage
  mns->rhsU   = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  mns->rhsV   = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  mns->rhsW   = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  mns->rhsP   = (dfloat*) calloc(Nlocal,sizeof(dfloat));
  mns->rhsPhi = (dfloat*) calloc(Nlocal,sizeof(dfloat));

  //additional field storage
  mns->NU   = (dfloat*) calloc(mns->NVfields*(mns->Nstages+1)*Ntotal,sizeof(dfloat));
  mns->LU   = (dfloat*) calloc(mns->NVfields*(mns->Nstages+1)*Ntotal,sizeof(dfloat));
  mns->GP   = (dfloat*) calloc(mns->NVfields*(mns->Nstages+1)*Ntotal,sizeof(dfloat));

  mns->GU   = (dfloat*) calloc(mns->NVfields*Ntotal*4,sizeof(dfloat));
  
  mns->rkU  = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
  mns->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  mns->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  mns->rkNU = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
  mns->rkLU = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
  mns->rkGP = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  mns->Vort = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
  mns->Div  = (dfloat*) calloc(              Nlocal,sizeof(dfloat));

  //extra storage for interpolated fields
  if(mns->elementType==HEXAHEDRA)
    mns->cU = (dfloat *) calloc(mns->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    mns->cU = mns->U;

  mns->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",mns->Nsubsteps);

  if(mns->Nsubsteps){
    mns->Ud    = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
    mns->Ue    = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
    mns->resU  = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));
    mns->rhsUd = (dfloat*) calloc(mns->NVfields*Ntotal,sizeof(dfloat));

    if(mns->elementType==HEXAHEDRA)
      mns->cUd = (dfloat *) calloc(mns->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      mns->cUd = mns->U;
  }

  dfloat rho  = 1.0 ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", mns->ubar);
  options.getArgs("VBAR", mns->vbar);
  if (mns->dim==3)
    options.getArgs("WBAR", mns->wbar);
  options.getArgs("PBAR", mns->pbar);
  options.getArgs("VISCOSITY", mns->nu);

  //Reynolds number
  mns->Re = mns->ubar/mns->nu;

  occa::kernelInfo kernelInfo;
  if(mns->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  occa::kernelInfo kernelInfoV  = kernelInfo;
  occa::kernelInfo kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo.addDefine("p_pbar", mns->pbar);
  kernelInfo.addDefine("p_ubar", mns->ubar);
  kernelInfo.addDefine("p_vbar", mns->vbar);
  kernelInfo.addDefine("p_wbar", mns->wbar);
  kernelInfo.addDefine("p_nu", mns->nu);

  kernelInfo.addDefine("p_NTfields", mns->NTfields);
  kernelInfo.addDefine("p_NVfields", mns->NVfields);
  kernelInfo.addDefine("p_NfacesNfp",  mesh->Nfaces*mesh->Nfp);
  kernelInfo.addDefine("p_Nstages",  mns->Nstages);
  kernelInfo.addDefine("p_SUBCYCLING",  mns->Nsubsteps);

  if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
    mns->ARKswitch = 1;   
  else 
    mns->ARKswitch = 0;

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo.addInclude((char*)boundaryHeaderFileName.c_str());

  mns->o_U   = mesh->device.malloc(mns->NVfields*mns->Nstages*Ntotal*sizeof(dfloat), mns->U);
  mns->o_P   = mesh->device.malloc(              mns->Nstages*Ntotal*sizeof(dfloat), mns->P);
  mns->o_Phi = mesh->device.malloc(                           Ntotal*sizeof(dfloat), mns->Phi);

  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);

  for (int r=0;r<size;r++) {
    if (r==rank) {
      if (mns->dim==2) 
        mns->setFlowFieldKernel =  mesh->device.buildKernelFromSource(DMNS "/okl/mnsSetFlowField2D.okl", "mnsSetFlowField2D", kernelInfo);  
      else
        mns->setFlowFieldKernel =  mesh->device.buildKernelFromSource(DMNS "/okl/mnsSetFlowField3D.okl", "mnsSetFlowField3D", kernelInfo);  
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  mns->startTime =0.0;
  options.getArgs("START TIME", mns->startTime);
  mns->setFlowFieldKernel(mesh->Nelements,
                          mns->startTime,
                          mesh->o_x,
                          mesh->o_y,
                          mesh->o_z,
                          mns->fieldOffset,
                          mns->o_U,
                          mns->o_P,
                          mns->o_Phi);
  mns->o_U.copyTo(mns->U);


  
  // mns->o_Phi.copyTo(mns->Phi);
  // char fname[BUFSIZ];
  // sprintf(fname,"TestPhi.vtu");
  // mnsPlotVTU(mns, fname);

 


  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      dfloat hest = 2./(sJ*invJ);

      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }
    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = mns->U[id+0*mns->fieldOffset];
      dfloat uyn = mns->U[id+1*mns->fieldOffset];
      dfloat uzn = 0.0;
      if (mns->dim==3) uzn = mns->U[id+2*mns->fieldOffset];


      //Squared maximum velocity
      dfloat numax;
      if (mns->dim==2)
        numax = uxn*uxn + uyn*uyn;
      else 
        numax = uxn*uxn + uyn*uyn + uzn*uzn;
      umax = mymax(umax, numax);
    }
  }

  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  
  dfloat cfl; 
  options.getArgs("CFL", cfl);
  dfloat dt     = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;

  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mns->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  options.getArgs("FINAL TIME", mns->finalTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    mns->NtimeSteps = mns->finalTime/mns->dt;

    if(mns->Nsubsteps){
      mns->dt         = mns->Nsubsteps*mns->dt;
      mns->NtimeSteps = mns->finalTime/mns->dt;
      mns->dt         = mns->finalTime/mns->NtimeSteps;
      mns->sdt        = mns->dt/mns->Nsubsteps;
    } else{
      mns->NtimeSteps = mns->finalTime/mns->dt;
      mns->dt         = mns->finalTime/mns->NtimeSteps;
    }
  }

  mns->dtMIN = 1E-2*mns->dt; //minumum allowed timestep

  if (rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  cfl);
    printf("dt = %g\n",   dt);
  }

  if (mns->Nsubsteps && rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", mns->dt, mns->sdt, mns->dt/mns->sdt);
  
  // Hold some inverses for kernels
  mns->inu = 1.0/mns->nu; 
  mns->idt = 1.0/mns->dt;
  
  mns->lambda = mns->g0 / (mns->dt * mns->nu);

  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", mns->outputStep);
  if (rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", mns->NtimeSteps,mns->outputStep, mns->dt);

#if 0
  //make option objects for elliptc solvers
  mns->vOptions = options;
  mns->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  mns->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  mns->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  mns->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  mns->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  mns->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  mns->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  mns->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  mns->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));

  mns->pOptions = options;
  mns->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  mns->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  mns->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  mns->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  mns->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  mns->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  mns->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  mns->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  mns->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));

  if (rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

  // SetUp Boundary Flags types for Elliptic Solve
  // bc = 1 -> wall
  // bc = 2 -> inflow
  // bc = 3 -> outflow
  // bc = 4 -> x-aligned slip
  // bc = 5 -> y-aligned slip
  // bc = 6 -> z-aligned slip
  int uBCType[7] = {0,1,1,2,1,2,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int vBCType[7] = {0,1,1,2,2,1,2}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int wBCType[7] = {0,1,1,2,2,2,1}; // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  int pBCType[7] = {0,2,2,1,2,2,2}; // bc=3 => outflow => Dirichlet => pBCType[3] = 1, etc.

  //Solver tolerances 
  mns->presTOL = 1E-8;
  mns->velTOL  = 1E-8;

  // Use third Order Velocity Solve: full rank should converge for low orders
  if (rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  mns->uSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mns->uSolver->mesh = mesh;
  mns->uSolver->options = mns->vOptions;
  mns->uSolver->dim = mns->dim;
  mns->uSolver->elementType = mns->elementType;
  mns->uSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mns->uSolver->BCType,uBCType,7*sizeof(int));
  ellipticSolveSetup(mns->uSolver, mns->lambda, kernelInfoV);

  mns->vSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mns->vSolver->mesh = mesh;
  mns->vSolver->options = mns->vOptions;
  mns->vSolver->dim = mns->dim;
  mns->vSolver->elementType = mns->elementType;
  mns->vSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mns->vSolver->BCType,vBCType,7*sizeof(int));
  ellipticSolveSetup(mns->vSolver, mns->lambda, kernelInfoV);

  if (mns->dim==3) {
    mns->wSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
    mns->wSolver->mesh = mesh;
    mns->wSolver->options = mns->vOptions;
    mns->wSolver->dim = mns->dim;
    mns->wSolver->elementType = mns->elementType;
    mns->wSolver->BCType = (int*) calloc(7,sizeof(int));
    memcpy(mns->wSolver->BCType,wBCType,7*sizeof(int));
    ellipticSolveSetup(mns->wSolver, mns->lambda, kernelInfoV);  
  }
  
  if (rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  mns->pSolver = (elliptic_t*) calloc(1, sizeof(elliptic_t));
  mns->pSolver->mesh = mesh;
  mns->pSolver->options = mns->pOptions;
  mns->pSolver->dim = mns->dim;
  mns->pSolver->elementType = mns->elementType;
  mns->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(mns->pSolver->BCType,pBCType,7*sizeof(int));
  ellipticSolveSetup(mns->pSolver, 0.0, kernelInfoP);


  //make node-wise boundary flags
  mns->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  mns->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) mns->VmapB[n+e*mesh->Np] = 1E9;
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
        for (int n=0;n<mesh->Nfp;n++) {
          int fid = mesh->faceNodes[n+f*mesh->Nfp];
          mns->VmapB[fid+e*mesh->Np] = mymin(bc,mns->VmapB[fid+e*mesh->Np]);
          mns->PmapB[fid+e*mesh->Np] = mymax(bc,mns->PmapB[fid+e*mesh->Np]);
        }
      }
    }
  }
  gsParallelGatherScatter(mesh->hostGsh, mns->VmapB, "int", "min"); 
  gsParallelGatherScatter(mesh->hostGsh, mns->PmapB, "int", "max"); 

  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {
    if (mns->VmapB[n] == 1E9) {
      mns->VmapB[n] = 0.;
    }
  }
  mns->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mns->VmapB);
  mns->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), mns->PmapB);

#endif

  kernelInfo.addDefine("p_blockSize", blockSize);
  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 256/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 256/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo.addDefine("p_maxNodesVolumeCub", maxNodesVolumeCub);
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo.addDefine("p_maxNodesSurfaceCub",maxNodesSurfaceCub);
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo.addDefine("p_cubNblockV",cubNblockV);
  kernelInfo.addDefine("p_cubNblockS",cubNblockS);

  // if (rank==0) {
  //   printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  //   printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

  //   printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  // }
  
  if (options.compareArgs("TIME INTEGRATOR", "ARK")) {
    mns->o_rkC  = mesh->device.malloc(         mns->Nrk*sizeof(dfloat),mns->rkC );
    mns->o_erkA = mesh->device.malloc(mns->Nrk*mns->Nrk*sizeof(dfloat),mns->erkA);
    mns->o_irkA = mesh->device.malloc(mns->Nrk*mns->Nrk*sizeof(dfloat),mns->irkA);
    mns->o_prkA = mesh->device.malloc(mns->Nrk*mns->Nrk*sizeof(dfloat),mns->prkA);
    mns->o_prkB = mesh->device.malloc(mns->Nrk*mns->Nrk*sizeof(dfloat),mns->prkB);
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    mns->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    mns->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    mns->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    mns->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    mns->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    mns->o_prkA = mns->o_extbdfC;
    mns->o_prkB = mns->o_extbdfC;
  }

  // MEMORY ALLOCATION
  mns->o_rhsU    = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsU);
  mns->o_rhsV    = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsV);
  mns->o_rhsW    = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsW);
  mns->o_rhsP    = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsP);
  mns->o_rhsPhi  = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsPhi);
  mns->o_resPhi  = mesh->device.malloc(Nlocal*sizeof(dfloat), mns->rhsPhi);

  mns->o_NU    = mesh->device.malloc(mns->NVfields*(mns->Nstages+1)*Ntotal*sizeof(dfloat), mns->NU);
  mns->o_LU    = mesh->device.malloc(mns->NVfields*(mns->Nstages+1)*Ntotal*sizeof(dfloat), mns->LU);
  mns->o_GP    = mesh->device.malloc(mns->NVfields*(mns->Nstages+1)*Ntotal*sizeof(dfloat), mns->GP);
  
  mns->o_GU    = mesh->device.malloc(mns->NVfields*Ntotal*4*sizeof(dfloat), mns->GU);
  
  mns->o_rkU   = mesh->device.malloc(mns->NVfields*Ntotal*sizeof(dfloat), mns->rkU);
  mns->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), mns->rkP);
  mns->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), mns->PI);
  
  mns->o_rkNU  = mesh->device.malloc(mns->NVfields*Ntotal*sizeof(dfloat), mns->rkNU);
  mns->o_rkLU  = mesh->device.malloc(mns->NVfields*Ntotal*sizeof(dfloat), mns->rkLU);
  mns->o_rkGP  = mesh->device.malloc(mns->NVfields*Ntotal*sizeof(dfloat), mns->rkGP);

  //storage for helmholtz solves
  mns->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  mns->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  mns->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  //plotting fields
  mns->o_Vort = mesh->device.malloc(mns->NVfields*Ntotal*sizeof(dfloat), mns->Vort);
  mns->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), mns->Div);

  if(mns->elementType==HEXAHEDRA)
    mns->o_cU = mesh->device.malloc(mns->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), mns->cU);
  else 
    mns->o_cU = mns->o_U;

  if(mesh->totalHaloPairs){//halo setup
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*(mns->NVfields)*sizeof(dfloat);
    occa::memory o_vsendBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    occa::memory o_vrecvBuffer = mesh->device.mappedAlloc(vHaloBytes, NULL);
    mns->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    mns->vSendBuffer = (dfloat*) o_vsendBuffer.getMappedPointer();
    mns->vRecvBuffer = (dfloat*) o_vrecvBuffer.getMappedPointer();

    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    mns->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);
    occa::memory o_psendBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    occa::memory o_precvBuffer = mesh->device.mappedAlloc(pHaloBytes, NULL);
    mns->pSendBuffer = (dfloat*) o_psendBuffer.getMappedPointer();
    mns->pRecvBuffer = (dfloat*) o_precvBuffer.getMappedPointer();

    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(mns->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat), NULL);
    mns->velocityHaloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer();
    mns->o_velocityHaloGatherTmp = mesh->device.malloc(mns->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat),  mns->velocityHaloGatherTmp);
  }

  // set kernel name suffix
  char *suffix;
  
  if(mns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mns->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mns->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<size;r++) {
    if (r==rank) {
      // Level Set kernels
      sprintf(fileName, DMNS "/okl/mnsLevelSet%s.okl", suffix);
      sprintf(kernelName, "mnsLevelSetCubatureVolume%s", suffix);
      mns->levelSetVolumeKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      sprintf(kernelName, "mnsLevelSetCubatureSurface%s",suffix);
      mns->levelSetSurfaceKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);
      sprintf(kernelName, "mnsLevelSetUpdate%s",suffix);
      mns->levelSetUpdateKernel = mesh->device.buildKernelFromSource(fileName,kernelName,kernelInfo);


      // sprintf(fileName, DINS "/okl/insHaloExchange.okl");
      // sprintf(kernelName, "insVelocityHaloExtract");
      // mns->velocityHaloExtractKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insVelocityHaloScatter");
      // mns->velocityHaloScatterKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insPressureHaloExtract");
      // ins->pressureHaloExtractKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insPressureHaloScatter");
      // ins->pressureHaloScatterKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================

      // sprintf(fileName, DINS "/okl/insAdvection%s.okl", suffix);
      // sprintf(kernelName, "insAdvectionCubatureVolume%s", suffix);
      // ins->advectionCubatureVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insAdvectionCubatureSurface%s", suffix);
      // ins->advectionCubatureSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insAdvectionVolume%s", suffix);
      // ins->advectionVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insAdvectionSurface%s", suffix);
      // ins->advectionSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insDiffusion%s.okl", suffix);
      // sprintf(kernelName, "insDiffusion%s", suffix);
      // ins->diffusionKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insDiffusionIpdg%s.okl", suffix);
      // sprintf(kernelName, "insDiffusionIpdg%s", suffix);
      // ins->diffusionIpdgKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insVelocityGradient%s.okl", suffix);
      // sprintf(kernelName, "insVelocityGradient%s", suffix);
      // ins->velocityGradientKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================

      // sprintf(fileName, DINS "/okl/insGradient%s.okl", suffix);
      // sprintf(kernelName, "insGradientVolume%s", suffix);
      // ins->gradientVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insGradientSurface%s", suffix);
      // ins->gradientSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insDivergence%s.okl", suffix);
      // sprintf(kernelName, "insDivergenceVolume%s", suffix);
      // ins->divergenceVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insDivergenceSurface%s", suffix);
      // ins->divergenceSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insVelocityRhs%s.okl", suffix);
      // if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
      //   sprintf(kernelName, "insVelocityRhsARK%s", suffix);
      // else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) 
      //   sprintf(kernelName, "insVelocityRhsEXTBDF%s", suffix);
      // ins->velocityRhsKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insVelocityBC%s.okl", suffix);
      // sprintf(kernelName, "insVelocityIpdgBC%s", suffix);
      // ins->velocityRhsIpdgBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insVelocityBC%s", suffix);
      // ins->velocityRhsBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insVelocityAddBC%s", suffix);
      // ins->velocityAddBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insPressureRhs%s.okl", suffix);
      // sprintf(kernelName, "insPressureRhs%s", suffix);
      // ins->pressureRhsKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insPressureBC%s.okl", suffix);
      // sprintf(kernelName, "insPressureIpdgBC%s", suffix);
      // ins->pressureRhsIpdgBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insPressureBC%s", suffix);
      // ins->pressureRhsBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(kernelName, "insPressureAddBC%s", suffix);
      // ins->pressureAddBCKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // // ===========================================================================

      // sprintf(fileName, DINS "/okl/insPressureUpdate.okl");
      // sprintf(kernelName, "insPressureUpdate");
      // ins->pressureUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insVelocityUpdate.okl");
      // sprintf(kernelName, "insVelocityUpdate");
      // ins->velocityUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);      

      // // ===========================================================================

      // sprintf(fileName, DINS "/okl/insVorticity%s.okl", suffix);
      // sprintf(kernelName, "insVorticity%s", suffix);
      // ins->vorticityKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);
    

      // if(ins->Nsubsteps){
      //   // Note that resU and resV can be replaced with already introduced buffer
      //   ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
      //   ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
      //   ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
      //   ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

      //   if(ins->elementType==HEXAHEDRA)
      //     ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
      //   else 
      //     ins->o_cUd = ins->o_Ud;

      //   sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
      //   sprintf(kernelName, "scaledAddwOffset");
      //   ins->scaledAddKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(fileName, DINS "/okl/insSubCycle%s.okl", suffix);
      //   sprintf(kernelName, "insSubCycleVolume%s", suffix);
      //   ins->subCycleVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(kernelName, "insSubCycleSurface%s", suffix);
      //   ins->subCycleSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(kernelName, "insSubCycleCubatureVolume%s", suffix);
      //   ins->subCycleCubatureVolumeKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(kernelName, "insSubCycleCubatureSurface%s", suffix);
      //   ins->subCycleCubatureSurfaceKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(fileName, DINS "/okl/insSubCycle.okl");
      //   sprintf(kernelName, "insSubCycleRKUpdate");
      //   ins->subCycleRKUpdateKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);

      //   sprintf(kernelName, "insSubCycleExt");
      //   ins->subCycleExtKernel =  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);
      // }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return mns;
}







