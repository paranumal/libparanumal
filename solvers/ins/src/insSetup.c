/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "ins.h"
#include "omp.h"
#include <unistd.h>

ins_t *insSetup(mesh_t *mesh, setupAide options){

  ins_t *ins = new ins_t();
  
  ins->mesh = mesh;
  ins->options = options;

  options.getArgs("MESH DIMENSION", ins->dim);
  options.getArgs("ELEMENT TYPE", ins->elementType);

  ins->NVfields = (ins->dim==3) ? 3:2; //  Total Number of Velocity Fields
  ins->NTfields = (ins->dim==3) ? 4:3; // Total Velocity + Pressure

  ins->frame = 0;
  ins->SNrk = 0;
  
  mesh->Nfields = 1; 

  ins->g0 =  1.0;

  if (options.compareArgs("TIME INTEGRATOR", "ARK1")) {
    ins->Nstages = 1;
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

    ins->Nrk = Nrk;
    ins->rkC = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->erkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->irkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkB = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));

    memcpy(ins->rkC, rkC, ins->Nrk*sizeof(dfloat));
    memcpy(ins->erkA, erkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkA, irkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkA, prkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkB, prkB, ins->Nrk*ins->Nrk*sizeof(dfloat));

    ins->g0 =  1.0;
    ins->embeddedRKFlag = 0; //no embedded method
  } else if (options.compareArgs("TIME INTEGRATOR", "ARK2")) {
    ins->Nstages = 2;
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

    ins->Nrk = Nrk;
    ins->rkC = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->erkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->irkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkB = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    

    memcpy(ins->rkC, rkC, ins->Nrk*sizeof(dfloat));
    memcpy(ins->erkA, erkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkA, irkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkA, prkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkB, prkB, ins->Nrk*ins->Nrk*sizeof(dfloat));
    
    ins->g0 =  1.0/gamma;
    ins->embeddedRKFlag = 0; //no embedded method
  } else if (options.compareArgs("TIME INTEGRATOR", "ARK3")) {
    ins->Nstages = 4;
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

    ins->Nrk = Nrk;
    ins->erkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->irkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkA = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->prkB = (dfloat*) calloc(ins->Nrk*ins->Nrk, sizeof(dfloat));
    ins->erkB = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->erkE = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->irkB = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->irkE = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));
    ins->rkC  = (dfloat*) calloc(ins->Nrk, sizeof(dfloat));


    memcpy(ins->erkA, erkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkA, irkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkA, prkA, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->prkB, prkB, ins->Nrk*ins->Nrk*sizeof(dfloat));
    memcpy(ins->erkB, erkB, ins->Nrk*sizeof(dfloat));
    memcpy(ins->erkE, erkE, ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkB, irkB, ins->Nrk*sizeof(dfloat));
    memcpy(ins->irkE, irkE, ins->Nrk*sizeof(dfloat));
    memcpy(ins->rkC, rkC, ins->Nrk*sizeof(dfloat));
    
    ins->g0 =  4055673282236.0/1767732205903.0;
    ins->embeddedRKFlag = 1; 
  } 


  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    ins->extbdfA = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfB = (dfloat*) calloc(3, sizeof(dfloat));
    ins->extbdfC = (dfloat*) calloc(3, sizeof(dfloat));

    ins->extC = (dfloat*) calloc(3, sizeof(dfloat));
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF1")) {
    ins->Nstages = 1;
    ins->temporalOrder = 1;
    ins->g0 = 1.0;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF2")) {
    ins->Nstages = 2;
    ins->temporalOrder = 2;
    ins->g0 = 1.5;
  } else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF3")) {
    ins->Nstages = 3;
    ins->temporalOrder = 3;
    ins->g0 = 11.f/6.f;
  }

  ins->readRestartFile = 0; 
  options.getArgs("RESTART FROM FILE", ins->readRestartFile);
  
  ins->writeRestartFile = 0; 
  options.getArgs("WRITE RESTART FILE", ins->writeRestartFile);



  dlong Nlocal = mesh->Np*mesh->Nelements;
  dlong Ntotal = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  
  ins->Ntotal = Ntotal;
  ins->fieldOffset = Ntotal;
  ins->Nblock = (Nlocal+blockSize-1)/blockSize;

  // compute samples of q at interpolation nodes
  ins->U     = (dfloat*) calloc(ins->NVfields*ins->Nstages*Ntotal,sizeof(dfloat));
  ins->P     = (dfloat*) calloc(              ins->Nstages*Ntotal,sizeof(dfloat));

  //rhs storage
  ins->rhsU  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsV  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsW  = (dfloat*) calloc(Ntotal,sizeof(dfloat));
  ins->rhsP  = (dfloat*) calloc(Ntotal,sizeof(dfloat));

  //additional field storage
  ins->NU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->LU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  ins->GP   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));

  ins->GU   = (dfloat*) calloc(ins->NVfields*Ntotal*4,sizeof(dfloat));
  
  ins->rkU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkP  = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  ins->PI   = (dfloat*) calloc(              Ntotal,sizeof(dfloat));
  
  ins->rkNU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkLU = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->rkGP = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

  //plotting fields
  ins->Vort = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
  ins->Div  = (dfloat*) calloc(              Nlocal,sizeof(dfloat));

  ins->FU   = (dfloat*) calloc(ins->NVfields*(ins->Nstages+1)*Ntotal,sizeof(dfloat));
  //extra storage for interpolated fields
  if(ins->elementType==HEXAHEDRA)
    ins->cU = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
  else 
    ins->cU = ins->U;

  ins->Nsubsteps = 0;
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    options.getArgs("SUBCYCLING STEPS",ins->Nsubsteps);

  if(ins->Nsubsteps){
    ins->Ud    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->Ue    = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->resU  = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));
    ins->rhsUd = (dfloat*) calloc(ins->NVfields*Ntotal,sizeof(dfloat));

    if(ins->elementType==HEXAHEDRA)
      ins->cUd = (dfloat *) calloc(ins->NVfields*mesh->Nelements*mesh->cubNp,sizeof(dfloat));
    else 
      ins->cUd = ins->U;

    // Prepare RK stages for Subcycling Part
    
    int Sorder = 4; // Defaulting to LSERK 4(5) 
    
    options.getArgs("SUBCYCLING TIME ORDER", Sorder);
    if(Sorder==2){
      ins->SNrk     = 2; 
      dfloat rka[2] = {0.0,     1.0 };
      dfloat rkb[2] = {0.5,     0.5 };
      dfloat rkc[2] = {0.0,     1.0 };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else if(Sorder ==3){
      // Using Williamson 3rd order scheme converted to low storage since the better truncation 
      ins->SNrk     = 3; 
      dfloat rka[3] = {0.0,     -5.0/9.0,  -153.0/128.0};
      dfloat rkb[3] = {1.0/3.0, 15.0/16.0,    8.0/15.0 };
      dfloat rkc[3] = {0.0,      1.0/3.0,     3.0/4.0  };
      //
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      //
      memcpy(ins->Srka, rka, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkb, rkb, ins->SNrk*sizeof(dfloat));
      memcpy(ins->Srkc, rkc, ins->SNrk*sizeof(dfloat));
    }else{
      ins->SNrk     = 5; 
      ins->Srka = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkb = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      ins->Srkc = (dfloat*) calloc(ins->SNrk, sizeof(dfloat));
      // Asumes initialized in mesh, can be moved here
      for(int rk=0; rk<ins->SNrk; rk++){
        ins->Srka[rk] = mesh->rka[rk]; 
        ins->Srkb[rk] = mesh->rkb[rk]; 
        ins->Srkc[rk] = mesh->rkc[rk]; 
      }
    }
  }

  dfloat rho  = 1.0 ;  // Give density for getting actual pressure in nondimensional solve
  dfloat g[3]; g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;  // No gravitational acceleration

  options.getArgs("UBAR", ins->ubar);
  options.getArgs("VBAR", ins->vbar);
  if (ins->dim==3)
    options.getArgs("WBAR", ins->wbar);
  options.getArgs("PBAR", ins->pbar);
  options.getArgs("VISCOSITY", ins->nu);

  //Reynolds number
  ins->Re = ins->ubar/ins->nu;

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(ins->dim==3){
    if(ins->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo); 
  } 
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  occa::properties kernelInfoV  = kernelInfo;
  occa::properties kernelInfoP  = kernelInfo;

  // ADD-DEFINES
  kernelInfo["defines/" "p_pbar"]= ins->pbar;
  kernelInfo["defines/" "p_ubar"]= ins->ubar;
  kernelInfo["defines/" "p_vbar"]= ins->vbar;
  kernelInfo["defines/" "p_wbar"]= ins->wbar;
  kernelInfo["defines/" "p_nu"]= ins->nu;

  kernelInfo["defines/" "p_NTfields"]= ins->NTfields;
  kernelInfo["defines/" "p_NVfields"]= ins->NVfields;
  kernelInfo["defines/" "p_NfacesNfp"]=  mesh->Nfaces*mesh->Nfp;
  kernelInfo["defines/" "p_Nstages"]=  ins->Nstages;
  kernelInfo["defines/" "p_SUBCYCLING"]=  ins->Nsubsteps;

  if(ins->elementType==QUADRILATERALS && mesh->dim==3){
    kernelInfo["defines/" "p_fainv"] = (dfloat) 0.0;
    kernelInfo["defines/" "p_invRadiusSq"] = (dfloat) 1./(mesh->sphereRadius*mesh->sphereRadius);
  }

  if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
    ins->ARKswitch = 1;   
  else 
    ins->ARKswitch = 0;

  //add boundary data to kernel info
  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  ins->o_U = mesh->device.malloc(ins->NVfields*ins->Nstages*Ntotal*sizeof(dfloat), ins->U);
  ins->o_P = mesh->device.malloc(              ins->Nstages*Ntotal*sizeof(dfloat), ins->P);

#if 0
  if (mesh->rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    occa::setVerboseCompilation(true);
  else 
    occa::setVerboseCompilation(false);
#endif
  

  if(options.compareArgs("INITIAL CONDITION", "BROWN-MINION") &&
     (ins->elementType == QUADRILATERALS && ins->dim==3)){
    printf("Setting up initial condition for BROWN-MINION test case...");
    insBrownMinionQuad3D(ins);
    ins->o_U.copyFrom(ins->U);
    ins->o_P.copyFrom(ins->P);
    printf("done\n");

    char fname[BUFSIZ];
    string outName;
    ins->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, ins->frame++);
    insPlotVTU(ins, fname);
  }else{

    for (int r=0;r<2;r++){
      if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {
	if (ins->dim==2) 
	  ins->setFlowFieldKernel =  mesh->device.buildKernel(DINS "/okl/insSetFlowField2D.okl", "insSetFlowField2D", kernelInfo);  
	else
	  ins->setFlowFieldKernel =  mesh->device.buildKernel(DINS "/okl/insSetFlowField3D.okl", "insSetFlowField3D", kernelInfo);  
      }
      MPI_Barrier(mesh->comm);
    }

    ins->startTime =0.0;
    options.getArgs("START TIME", ins->startTime);
    ins->setFlowFieldKernel(mesh->Nelements,
			    ins->startTime,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    ins->fieldOffset,
			    ins->o_U,
			    ins->o_P);
    ins->o_U.copyTo(ins->U);
    
  }
  
  // set time step
  dfloat hmin = 1e9, hmax = 0;
  dfloat umax = 0;
  for(dlong e=0;e<mesh->Nelements;++e){

    if(ins->elementType==TRIANGLES || ins->elementType == TETRAHEDRA){
      for(int f=0;f<mesh->Nfaces;++f){
	dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
	dfloat sJ   = mesh->sgeo[sid + SJID];
	dfloat invJ = mesh->sgeo[sid + IJID];

	dfloat hest = 2./(sJ*invJ);

	hmin = mymin(hmin, hest);
	hmax = mymax(hmax, hest);
      }
    }else{
      for(int f=0;f<mesh->Nfaces;++f){
	for(int n=0; n<mesh->Nfp; n++){
	  dlong sid = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n);
	  dfloat sJ   = mesh->sgeo[sid + SJID];
	  dfloat invJ = mesh->sgeo[sid + IJID];

	  dfloat hest = 2./(sJ*invJ);

	  hmin = mymin(hmin, hest);
	  hmax = mymax(hmax, hest);
	}
      }
    }

    // dfloat maxMagVecLoc = 0;

    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;
      dfloat t = 0;
      dfloat uxn = ins->U[id+0*ins->fieldOffset];
      dfloat uyn = ins->U[id+1*ins->fieldOffset];
      dfloat uzn = 0.0;
      if (ins->dim==3) uzn = ins->U[id+2*ins->fieldOffset];


      //Squared maximum velocity
      dfloat numax;
      if (ins->dim==2)
	numax = uxn*uxn + uyn*uyn;
      else 
	numax = uxn*uxn + uyn*uyn + uzn*uzn;

      umax = mymax(umax, numax);
    }
  }

  // Maximum Velocity
  umax = sqrt(umax);
  dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
  printf("magVel = %lf\n", magVel);
  
  options.getArgs("CFL", ins->cfl);
  dfloat dt     = ins->cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
  
  ins->dtAdaptStep = 0; 
  options.getArgs("TSTEPS FOR TIME STEP ADAPT", ins->dtAdaptStep);

  // over ride if DT is in the file
  if(options.getArgs("DT", dt)){
    if(mesh->rank==0){
      printf("NOTE: USING DT %lg FROM OPTIONS FILE\n", dt);
    }
  }
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dti), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  // save initial time-step estimate 
  ins->dt = ins->dti; 

  options.getArgs("FINAL TIME", ins->finalTime);
  options.getArgs("START TIME", ins->startTime);
  
  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;

    if(ins->Nsubsteps){
      ins->dt         = ins->Nsubsteps*ins->dt;
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
      ins->sdt        = ins->dt/ins->Nsubsteps;
    } else{
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
    }
  }

  ins->dtMIN = 1E-2*ins->dt; //minumum allowed timestep

  if (mesh->rank==0) {
    printf("hmin = %g\n", hmin);
    printf("hmax = %g\n", hmax);
    printf("cfl = %g\n",  ins->cfl);
    printf("dt = %g\n",   dt);
  }

  if (ins->Nsubsteps && mesh->rank==0) printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);
  
  // Hold some inverses for kernels
  ins->inu = 1.0/ins->nu; 
  ins->idt = 1.0/ins->dt;
  
  ins->lambda = ins->g0 / (ins->dt * ins->nu);

  options.getArgs("TSTEPS FOR SOLUTION OUTPUT", ins->outputStep);
  if (mesh->rank==0) printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->outputStep, ins->dt);

  ins->outputForceStep = 0;
  options.getArgs("TSTEPS FOR FORCE OUTPUT", ins->outputForceStep);

  ins->isoNfields  = 1;   //1 + (ins->dim) + (1 + ins->dim) ; // p, u.v,w, vort_x, vort_y, vort_z, wort_mag 
  ins->isoMaxNtris = 1.E7; 
  
  //!!!!! Isosurface Setup !! remove those to seoperate library later !!!!
  if(ins->dim==3 && options.compareArgs("OUTPUT TYPE", "ISO"))
    {
    
      // Only one field is exported for iso-surface to reduce the file size

      //
      options.getArgs("ISOSURFACE FIELD ID", ins->isoField); 
      options.getArgs("ISOSURFACE COLOR ID", ins->isoColorField); 
      options.getArgs("ISOSURFACE LEVEL NUMBER", ins->isoNlevels);
      options.getArgs("ISOSURFACE CONTOUR MAX", ins->isoMaxVal); 
      options.getArgs("ISOSURFACE CONTOUR MIN", ins->isoMinVal);

      ins->isoMax    = (ins->dim + ins->isoNfields)*3*ins->isoMaxNtris;
      ins->isoNtris  = (int*) calloc(1, sizeof(int));
      ins->isoq      = (dfloat*) calloc(ins->isoMax, sizeof(dfloat)); 

      ins->o_isoq      = mesh->device.malloc(ins->isoMax*sizeof(dfloat), ins->isoq);
      ins->o_isoNtris  = mesh->device.malloc(1*sizeof(int), ins->isoNtris);

      // Create all contour levels
      dfloat *isoLevels = (dfloat*) calloc(ins->isoNlevels, sizeof(dfloat));
      for(int l=0;l<ins->isoNlevels;++l)
	isoLevels[l] = ins->isoMinVal + (ins->isoMaxVal-ins->isoMinVal)*l/(dfloat)(ins->isoNlevels-1);



      // GROUP LEVELS of ISOCONTOURS

      int levelsInGroup = 0; 
      options.getArgs("ISOSURFACE GROUP NUMBER", levelsInGroup);

      if(levelsInGroup==0) {printf("Number of levels in each group can not be zero!!!\n");  exit(EXIT_FAILURE);} 
      if(levelsInGroup){

	// Number of groups for isosurfaces
	ins->isoGNgroups        = ins->isoNlevels/(levelsInGroup);  
	if(ins->isoNlevels%(levelsInGroup))
	  ins->isoGNgroups++; 

	ins->isoGNlevels        = (int *) calloc(ins->isoGNgroups,sizeof(int));
	ins->isoGLvalues        = (dfloat **) calloc(ins->isoGNgroups,sizeof(dfloat*));

	for(int gr =0; gr<ins->isoGNgroups; gr++)
	  {
	    int nlevels = (gr+1)*levelsInGroup > ins->isoNlevels ? (ins->isoNlevels%levelsInGroup) : levelsInGroup;  
	    ins->isoGNlevels[gr] = nlevels;  
	    printf("Isosurface Group %d has %d levels\n", gr, ins->isoGNlevels[gr]);
	  }

	// Allocate memory for levels in each group
	for (int gr =0;gr<ins->isoGNgroups;gr++)
	  ins->isoGLvalues[gr] = (dfloat *) calloc(ins->isoGNlevels[gr],sizeof(dfloat));

	int sk = 0; 
	for (int gr =0;gr<ins->isoGNgroups;gr++){
	  printf("Isosurface Group %d Values\n", gr);        
	  for (int l=0;l<ins->isoGNlevels[gr];l++){
	    ins->isoGLvalues[gr][l] = isoLevels[sk + l];
	    printf("%.4f\t", ins->isoGLvalues[gr][l]);
	  }
	  printf("\n");
	  sk += ins->isoGNlevels[gr]; 
	}

	// Create levels for each group
	ins->o_isoGLvalues     = (occa::memory *) malloc(ins->isoGNgroups*sizeof(occa::memory));
	for (int gr =0;gr<ins->isoGNgroups;gr++)
	  ins->o_isoGLvalues[gr] = mesh->device.malloc(ins->isoGNlevels[gr]*sizeof(dfloat),ins->isoGLvalues[gr]);
    
      }

      // Interpolation operators form Np to PlotNp (equisapaced nodes of order >N generally)
      dfloat *plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
      for(int n=0;n<mesh->plotNp;++n){
	for(int m=0;m<mesh->Np;++m){
	  plotInterp[n+m*mesh->plotNp] = mesh->plotInterp[n*mesh->Np+m];
	}
      }
      ins->o_plotInterp = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(dfloat), plotInterp);

      // EToV for local triangulation
      int *plotEToV = (int*) calloc(mesh->plotNp*mesh->Np, sizeof(int));
      for(int n=0;n<mesh->plotNelements;++n){
	for(int m=0;m<mesh->plotNverts;++m){
	  plotEToV[n+m*mesh->plotNelements] = mesh->plotEToV[n*mesh->plotNverts+m];
	}
      }
      ins->o_plotEToV = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(int), plotEToV);


    }

  //make option objects for elliptc solvers
  ins->vOptions = options;
  ins->vOptions.setArgs("KRYLOV SOLVER",        options.getArgs("VELOCITY KRYLOV SOLVER"));
  ins->vOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("VELOCITY SOLVER TOLERANCE"));
  ins->vOptions.setArgs("DISCRETIZATION",       options.getArgs("VELOCITY DISCRETIZATION"));
  ins->vOptions.setArgs("BASIS",                options.getArgs("VELOCITY BASIS"));
  ins->vOptions.setArgs("PRECONDITIONER",       options.getArgs("VELOCITY PRECONDITIONER"));
  ins->vOptions.setArgs("MULTIGRID COARSENING", options.getArgs("VELOCITY MULTIGRID COARSENING"));
  ins->vOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("VELOCITY MULTIGRID SMOOTHER"));
  ins->vOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("VELOCITY MULTIGRID CHEBYSHEV DEGREE"));

  ins->vOptions.setArgs("PARALMOND CYCLE",      options.getArgs("VELOCITY PARALMOND CYCLE"));
  ins->vOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("VELOCITY PARALMOND SMOOTHER"));
  ins->vOptions.setArgs("PARALMOND PARTITION",  options.getArgs("VELOCITY PARALMOND PARTITION"));
  ins->vOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("VELOCITY PARALMOND CHEBYSHEV DEGREE"));

  ins->vOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("VELOCITY PARALMOND AGGREGATION STRATEGY"));
  
  ins->vOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->vOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  ///  std::cout << "vOptions: " << ins->vOptions << std::endl;
  
  ins->pOptions = options;
  ins->pOptions.setArgs("KRYLOV SOLVER",        options.getArgs("PRESSURE KRYLOV SOLVER"));
  ins->pOptions.setArgs("SOLVER TOLERANCE",     options.getArgs("PRESSURE SOLVER TOLERANCE"));
  ins->pOptions.setArgs("DISCRETIZATION",       options.getArgs("PRESSURE DISCRETIZATION"));
  ins->pOptions.setArgs("BASIS",                options.getArgs("PRESSURE BASIS"));
  ins->pOptions.setArgs("PRECONDITIONER",       options.getArgs("PRESSURE PRECONDITIONER"));
  ins->pOptions.setArgs("MULTIGRID COARSENING", options.getArgs("PRESSURE MULTIGRID COARSENING"));
  ins->pOptions.setArgs("MULTIGRID SMOOTHER",   options.getArgs("PRESSURE MULTIGRID SMOOTHER"));
  ins->pOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE",  options.getArgs("PRESSURE MULTIGRID CHEBYSHEV DEGREE"));

  ins->pOptions.setArgs("PARALMOND CYCLE",      options.getArgs("PRESSURE PARALMOND CYCLE"));
  ins->pOptions.setArgs("PARALMOND SMOOTHER",   options.getArgs("PRESSURE PARALMOND SMOOTHER"));
  ins->pOptions.setArgs("PARALMOND PARTITION",  options.getArgs("PRESSURE PARALMOND PARTITION"));
  ins->pOptions.setArgs("PARALMOND CHEBYSHEV DEGREE",  options.getArgs("PRESSURE PARALMOND CHEBYSHEV DEGREE"));

  ins->pOptions.setArgs("PARALMOND AGGREGATION STRATEGY", options.getArgs("PRESSURE PARALMOND AGGREGATION STRATEGY"));
  
  ins->pOptions.setArgs("DEBUG ENABLE OGS", "1");
  ins->pOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");

  ///  std::cout << "pOptions: " << ins->pOptions << std::endl;  
  
  if (mesh->rank==0) printf("==================ELLIPTIC SOLVE SETUP=========================\n");

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
  ins->presTOL = 1E-8;
  ins->velTOL  = 1E-8;


  // Use third Order Velocity Solve: full rank should converge for low orders
  if (mesh->rank==0) printf("==================VELOCITY SOLVE SETUP=========================\n");

  ins->uSolver = new elliptic_t();
  ins->uSolver->mesh = mesh;
  ins->uSolver->options = ins->vOptions;
  ins->uSolver->dim = ins->dim;
  ins->uSolver->elementType = ins->elementType;
  ins->uSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->uSolver->BCType,uBCType,7*sizeof(int));

  ellipticSolveSetup(ins->uSolver, ins->lambda, kernelInfoV); 

  ins->vSolver = new elliptic_t();
  ins->vSolver->mesh = mesh;
  ins->vSolver->options = ins->vOptions;
  ins->vSolver->dim = ins->dim;
  ins->vSolver->elementType = ins->elementType;
  ins->vSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->vSolver->BCType,vBCType,7*sizeof(int));
  
  ellipticSolveSetup(ins->vSolver, ins->lambda, kernelInfoV); //!!!!!

  
  if (ins->dim==3) {
    ins->wSolver = new elliptic_t();
    ins->wSolver->mesh = mesh;
    ins->wSolver->options = ins->vOptions;
    ins->wSolver->dim = ins->dim;
    ins->wSolver->elementType = ins->elementType;
    ins->wSolver->BCType = (int*) calloc(7,sizeof(int));
    memcpy(ins->wSolver->BCType,wBCType,7*sizeof(int));
    
    ellipticSolveSetup(ins->wSolver, ins->lambda, kernelInfoV);  //!!!!!
  }
  
  if (mesh->rank==0) printf("==================PRESSURE SOLVE SETUP=========================\n");
  ins->pSolver = new elliptic_t();
  ins->pSolver->mesh = mesh;
  ins->pSolver->options = ins->pOptions;
  ins->pSolver->dim = ins->dim;
  ins->pSolver->elementType = ins->elementType;
  ins->pSolver->BCType = (int*) calloc(7,sizeof(int));
  memcpy(ins->pSolver->BCType,pBCType,7*sizeof(int));
  
  ellipticSolveSetup(ins->pSolver, 0.0, kernelInfoP); //!!!!

  // TW: this code needs to be re-evaluated from here ====>
  //make node-wise boundary flags
  dfloat largeNumber = 1<<20;
  ins->VmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  ins->PmapB = (int *) calloc(mesh->Nelements*mesh->Np,sizeof(int));
  for (int e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) ins->VmapB[n+e*mesh->Np] = largeNumber;
  }

  for (int e=0;e<mesh->Nelements;e++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      int bc = mesh->EToB[f+e*mesh->Nfaces];
      if (bc>0) {
	for (int n=0;n<mesh->Nfp;n++) {
	  int fid = mesh->faceNodes[n+f*mesh->Nfp];
	  ins->VmapB[fid+e*mesh->Np] = mymin(bc,ins->VmapB[fid+e*mesh->Np]);
	  ins->PmapB[fid+e*mesh->Np] = mymax(bc,ins->PmapB[fid+e*mesh->Np]);
	}
      }
    }
  }

  // this is potentially not good for Neumann (since it will )
  ogsGatherScatter(ins->VmapB, ogsInt, ogsMin, mesh->ogs); // !!!!!!!!!!!!!!!!!
  ogsGatherScatter(ins->PmapB, ogsInt, ogsMax, mesh->ogs); // !!!!!!!!!!!!!!!!!
#if 1
  for (int n=0;n<mesh->Nelements*mesh->Np;n++) {

    if (ins->VmapB[n] == largeNumber) {
      ins->VmapB[n] = 0.;
    }
  }
#endif
 
  ins->o_VmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->VmapB);
  ins->o_PmapB = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(int), ins->PmapB);

  // TW: this code needs to be re-evaluated to here <======
  // and the kernels that use VmapB, PmapB
  

  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  if(options.compareArgs("TIME INTEGRATOR", "EXTBDF"))
    kernelInfo["defines/" "p_EXTBDF"]= 1;
  else
    kernelInfo["defines/" "p_EXTBDF"]= 0;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = mymax(1,256/mesh->Np); // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = mymax(1,256/maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;


  int maxNodesVolumeCub = mymax(mesh->cubNp,mesh->Np);  
  kernelInfo["defines/" "p_maxNodesVolumeCub"]= maxNodesVolumeCub;
  int cubNblockV = mymax(1,256/maxNodesVolumeCub);
  //
  int maxNodesSurfaceCub = mymax(mesh->Np, mymax(mesh->Nfaces*mesh->Nfp, mesh->Nfaces*mesh->intNfp));
  kernelInfo["defines/" "p_maxNodesSurfaceCub"]=maxNodesSurfaceCub;
  int cubNblockS = mymax(256/maxNodesSurfaceCub,1); // works for CUDA
  //
  kernelInfo["defines/" "p_cubNblockV"]=cubNblockV;
  kernelInfo["defines/" "p_cubNblockS"]=cubNblockS;

  
  // IsoSurface related
  if(ins->dim==3 && ins->elementType != QUADRILATERALS ){
    kernelInfo["defines/" "p_isoNfields"]= ins->isoNfields;
    // Define Isosurface Area Tolerance
    kernelInfo["defines/" "p_triAreaTol"]= (dfloat) 1.0E-16;

    kernelInfo["defines/" "p_dim"]= ins->dim;
    kernelInfo["defines/" "p_plotNp"]= mesh->plotNp;
    kernelInfo["defines/" "p_plotNelements"]= mesh->plotNelements;
    
    int plotNthreads = mymax(mesh->Np, mymax(mesh->plotNp, mesh->plotNelements));
    kernelInfo["defines/" "p_plotNthreads"]= plotNthreads;
  } 





  // if (mesh->rank==0) {
  //   printf("maxNodes: %d \t  NblockV: %d \t NblockS: %d  \n", maxNodes, NblockV, NblockS);
  //   printf("maxNodesVolCub: %d \t maxNodesSurCub: %d \t NblockVCub: %d \t NblockSCub: %d  \n", maxNodesVolumeCub,maxNodesSurfaceCub, cubNblockV, cubNblockS);

  //   printf("Np: %d \t Ncub: %d \n", mesh->Np, mesh->cubNp);
  // }
  
  if (options.compareArgs("TIME INTEGRATOR", "ARK")) {
    ins->o_rkC  = mesh->device.malloc(         ins->Nrk*sizeof(dfloat),ins->rkC );
    ins->o_erkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->erkA);
    ins->o_irkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->irkA);
    ins->o_prkA = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->prkA);
    ins->o_prkB = mesh->device.malloc(ins->Nrk*ins->Nrk*sizeof(dfloat),ins->prkB);
  }

  if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    dfloat rkC[4] = {1.0, 0.0, -1.0, -2.0};

    ins->o_rkC  = mesh->device.malloc(4*sizeof(dfloat),rkC);
    ins->o_extbdfA = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfB = mesh->device.malloc(3*sizeof(dfloat));
    ins->o_extbdfC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_extC = mesh->device.malloc(3*sizeof(dfloat)); 

    ins->o_prkA = ins->o_extbdfC;
    ins->o_prkB = ins->o_extbdfC;
  }

  // MEMORY ALLOCATION
  ins->o_rhsU  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsU);
  ins->o_rhsV  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsV);
  ins->o_rhsW  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsW);
  ins->o_rhsP  = mesh->device.malloc(Ntotal*sizeof(dfloat), ins->rhsP);

  ins->o_NU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->NU);
  ins->o_LU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->LU);
  ins->o_GP    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->GP);
  
  ins->o_GU    = mesh->device.malloc(ins->NVfields*Ntotal*4*sizeof(dfloat), ins->GU);
  
  ins->o_rkU   = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkU);
  ins->o_rkP   = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->rkP);
  ins->o_PI    = mesh->device.malloc(              Ntotal*sizeof(dfloat), ins->PI);
  
  ins->o_rkNU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkNU);
  ins->o_rkLU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkLU);
  ins->o_rkGP  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rkGP);

  //storage for helmholtz solves
  ins->o_UH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_VH = mesh->device.malloc(Ntotal*sizeof(dfloat));
  ins->o_WH = mesh->device.malloc(Ntotal*sizeof(dfloat));

  //plotting fields
  ins->o_Vort = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Vort);
  ins->o_Div  = mesh->device.malloc(              Nlocal*sizeof(dfloat), ins->Div);

  ins->o_FU    = mesh->device.malloc(ins->NVfields*(ins->Nstages+1)*Ntotal*sizeof(dfloat), ins->FU);

  if(ins->elementType==HEXAHEDRA) // !!!! check that
    ins->o_cU = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cU);
  else 
    ins->o_cU = ins->o_U;

  if(mesh->totalHaloPairs){//halo setup
    //    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
    dlong vHaloBytes = mesh->totalHaloPairs*mesh->Np*ins->NVfields*sizeof(dfloat);
    dlong pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    dlong vGatherBytes = ins->NVfields*mesh->ogs->NhaloGather*sizeof(dfloat);

    ins->o_vHaloBuffer = mesh->device.malloc(vHaloBytes);
    ins->o_pHaloBuffer = mesh->device.malloc(pHaloBytes);

    //    occa::memory o_vSendBuffer, o_vRecvBuffer, o_pSendBuffer, o_pRecvBuffer, o_gatherTmpPinned;
#if 1
    ins->vSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vSendBuffer, ins->h_vSendBuffer);
    ins->vRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, vHaloBytes, NULL, ins->o_vRecvBuffer, ins->h_vRecvBuffer);

    ins->pSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pSendBuffer, ins->h_pSendBuffer);
    ins->pRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, pHaloBytes, NULL, ins->o_pRecvBuffer, ins->h_pRecvBuffer);

    ins->velocityHaloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, vGatherBytes, NULL, ins->o_gatherTmpPinned, ins->h_gatherTmpPinned);
    ins->o_velocityHaloGatherTmp = mesh->device.malloc(vGatherBytes,  ins->velocityHaloGatherTmp);
#else
    ins->vSendBuffer = (dfloat*) malloc(vHaloBytes);
    ins->vRecvBuffer = (dfloat*) malloc(vHaloBytes);
    ins->pSendBuffer = (dfloat*) malloc(pHaloBytes);
    ins->pRecvBuffer = (dfloat*) malloc(pHaloBytes);
#endif
  }

  // set kernel name suffix
  char *suffix;
  
  if(ins->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(ins->elementType==QUADRILATERALS){
    if(ins->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D"); 
  }
  if(ins->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(ins->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  
  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      
      //      std::cout << "NEW " << " " <<  kernelInfo << std::endl;
      
      sprintf(fileName, DINS "/okl/insHaloExchange.okl");
      sprintf(kernelName, "insVelocityHaloExtract");
      ins->velocityHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insVelocityHaloScatter");
      ins->velocityHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloExtract");
      ins->pressureHaloExtractKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insPressureHaloScatter");
      ins->pressureHaloScatterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // --
      if(ins->dim==3 && ins->elementType==QUADRILATERALS){
      	sprintf(fileName, DINS "/okl/insConstrainQuad3D.okl");
      	sprintf(kernelName, "insConstrainQuad3D");
      	ins->constrainKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }
      
      // ===========================================================================

      sprintf(fileName, DINS "/okl/insAdvection%s.okl", suffix);

      // needed to be implemented
      sprintf(kernelName, "insAdvectionCubatureVolume%s", suffix);
      ins->advectionCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionCubatureSurface%s", suffix);
      ins->advectionCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionVolume%s", suffix);
      ins->advectionVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insAdvectionSurface%s", suffix);
      ins->advectionSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insStrongAdvectionVolume%s", suffix);
      ins->advectionStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      // ===========================================================================
      if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION")){
        sprintf(fileName, DINS "/okl/insFilter.okl");

        sprintf(kernelName, "insFilterRT%s", suffix);
        ins->filterKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }


      // ===========================================================================
      
      // sprintf(fileName, DINS "/okl/insDiffusion%s.okl", suffix);
      // sprintf(kernelName, "insDiffusion%s", suffix);
      // ins->diffusionKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insDiffusionIpdg%s.okl", suffix);
      // sprintf(kernelName, "insDiffusionIpdg%s", suffix);
      // ins->diffusionIpdgKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // sprintf(fileName, DINS "/okl/insVelocityGradient%s.okl", suffix);
      // sprintf(kernelName, "insVelocityGradient%s", suffix);
      // ins->velocityGradientKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // // ===========================================================================

      sprintf(fileName, DINS "/okl/insGradient%s.okl", suffix);
      sprintf(kernelName, "insGradientVolume%s", suffix);
      ins->gradientVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(kernelName, "insGradientSurface%s", suffix);
      ins->gradientSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insDivergence%s.okl", suffix);

      sprintf(kernelName, "insDivergenceVolume%s", suffix);
      ins->divergenceVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      sprintf(kernelName, "insDivergenceSurface%s", suffix);
      ins->divergenceSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      
      // AK. Note that currently only HEX and QUAD is implemented in weak form !!!!!!!
      if(ins->elementType==HEXAHEDRA || (ins->elementType==QUADRILATERALS && ins->dim==2)){
        sprintf(kernelName, "insStrongDivergenceVolume%s", suffix);
        ins->divergenceStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }else{ // implement other divergence operators in weak form also!!!!
        sprintf(kernelName, "insDivergenceVolume%s", suffix);
        ins->divergenceStrongVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo); 
      }

      
      // ===========================================================================
      
      sprintf(fileName, DINS "/okl/insVelocityRhs%s.okl", suffix);
      if (options.compareArgs("TIME INTEGRATOR", "ARK")) 
	sprintf(kernelName, "insVelocityRhsARK%s", suffix); 
      else if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) 
	sprintf(kernelName, "insVelocityRhsEXTBDF%s", suffix);
      ins->velocityRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);



      if(!(ins->dim==3 && ins->elementType==QUADRILATERALS) ){
	sprintf(fileName, DINS "/okl/insVelocityBC%s.okl", suffix);
	sprintf(kernelName, "insVelocityIpdgBC%s", suffix);
	ins->velocityRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insVelocityBC%s", suffix);
	ins->velocityRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insVelocityAddBC%s", suffix);
	ins->velocityAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // ===========================================================================
      
      // Dont forget to modify!!!!!!!!
      sprintf(fileName, DINS "/okl/insPressureRhs%s.okl", suffix);
      sprintf(kernelName, "insPressureRhs%s", suffix);
      ins->pressureRhsKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      
      if(!(ins->dim==3 && ins->elementType==QUADRILATERALS) ){
	sprintf(fileName, DINS "/okl/insPressureBC%s.okl", suffix);
	sprintf(kernelName, "insPressureIpdgBC%s", suffix);
	ins->pressureRhsIpdgBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insPressureBC%s", suffix);
	ins->pressureRhsBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insPressureAddBC%s", suffix);
	ins->pressureAddBCKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insPressureUpdate.okl");
      sprintf(kernelName, "insPressureUpdate");
      ins->pressureUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      sprintf(fileName, DINS "/okl/insVelocityUpdate.okl");
      sprintf(kernelName, "insVelocityUpdate");
      ins->velocityUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);      

      // ===========================================================================

      sprintf(fileName, DINS "/okl/insVorticity%s.okl", suffix);
      sprintf(kernelName, "insVorticity%s", suffix);
      ins->vorticityKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
    
      // ===========================================================================
      if(ins->dim==3 && ins->options.compareArgs("OUTPUT TYPE","ISO")){
	sprintf(fileName, DINS "/okl/insIsoSurface3D.okl");
	sprintf(kernelName, "insIsoSurface3D");

	ins->isoSurfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);  
      }
      

      // Not implemented for Quad 3D yet !!!!!!!!!!
      if(ins->Nsubsteps){
	// Note that resU and resV can be replaced with already introduced buffer
	ins->o_Ue    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ue);
	ins->o_Ud    = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->Ud);
	ins->o_resU  = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->resU);
	ins->o_rhsUd = mesh->device.malloc(ins->NVfields*Ntotal*sizeof(dfloat), ins->rhsUd);

	if(ins->elementType==HEXAHEDRA)
	  ins->o_cUd = mesh->device.malloc(ins->NVfields*mesh->Nelements*mesh->cubNp*sizeof(dfloat), ins->cUd);
	else 
	  ins->o_cUd = ins->o_Ud;
      
	sprintf(fileName, DHOLMES "/okl/scaledAdd.okl");
	sprintf(kernelName, "scaledAddwOffset");
	ins->scaledAddKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(fileName, DINS "/okl/insSubCycle%s.okl", suffix);
	sprintf(kernelName, "insSubCycleVolume%s", suffix);
	ins->subCycleVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insSubCycleSurface%s", suffix);
	ins->subCycleSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insSubCycleCubatureVolume%s", suffix);
	ins->subCycleCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  
  if(ins->elementType==HEXAHEDRA){
	sprintf(kernelName, "insSubCycleNekCubatureVolume%s", suffix);
	ins->subCycleNekCubatureVolumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insSubCycleCubatureSurface%s", suffix);
	ins->subCycleCubatureSurfaceKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
  }

	sprintf(fileName, DINS "/okl/insSubCycle.okl");
	sprintf(kernelName, "insSubCycleRKUpdate");
	ins->subCycleRKUpdateKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

	sprintf(kernelName, "insSubCycleExt");
	ins->subCycleExtKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);
      }

      ins->haloGetKernel =
	mesh->device.buildKernel(DINS "/okl/insHalo.okl", "insHaloGet", kernelInfo);
      ins->haloPutKernel =
	mesh->device.buildKernel(DINS "/okl/insHalo.okl", "insHaloPut", kernelInfo);


    sprintf(fileName, DHOLMES "/okl/addScalar.okl");
    sprintf(kernelName, "setScalar");
    ins->setScalarKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

      
    }
    MPI_Barrier(mesh->comm);
  }

  // build lumped mass matrix for NEK
  dfloat *lumpedMassMatrix =
    (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
  dfloat *copyLumpedMassMatrix =
    (dfloat*) calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));

  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      lumpedMassMatrix[e*mesh->Np+n] =
	mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
      copyLumpedMassMatrix[e*mesh->Np+n] =
	mesh->vgeo[e*mesh->Np*mesh->Nvgeo+JWID*mesh->Np+n];
    }
  }
  
  ogsGatherScatter(lumpedMassMatrix, ogsDfloat, ogsAdd, ins->vSolver->ogs); // !!!!!!!!!!!!!!!!!

  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    //    lumpedMassMatrix[n] = copyLumpedMassMatrix[n]/lumpedMassMatrix[n];
    lumpedMassMatrix[n] = 1./lumpedMassMatrix[n];
  }

  ins->o_invLumpedMassMatrix = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),
						   lumpedMassMatrix);

 
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
  insFilterSetup(ins); 
  
  return ins;
}







