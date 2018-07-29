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

#include "elliptic.h"

typedef struct{

  dlong localId;
  hlong baseId;
  int haloFlag;

} preconGatherInfo_t;

int parallelCompareBaseId(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;
}

void ellipticSetupSmootherLocalPatch(elliptic_t *elliptic, precon_t *precon, 
                                      agmgLevel *level, dfloat lambda, 
                                      dfloat rateTolerance) {

  dfloat *invAP;
  dlong Npatches;
  dlong *patchesIndex;

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int NpP = mesh->Np;

  //initialize the full inverse operators on each 4 element patch
  ellipticBuildLocalPatches(elliptic, lambda, rateTolerance, &Npatches, &patchesIndex, &invAP);

  precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
  precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(dlong), patchesIndex);

  dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
  for (dlong e=0;e<mesh->Nelements;e++) {
    invDegree[e] = 1.0;
  }
  precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);

  level->device_smoother = LocalPatch;

  //estimate the max eigenvalue of S*A
  dfloat rho = maxEigSmoothAx(elliptic, level);

  if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {

    level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

    level->smoother_params[0] = rho;
    level->smoother_params[1] = rho/10.;

  } else {

    //set the stabilty weight (jacobi-type interation)
    dfloat weight = (4./3.)/rho;

    for (dlong e=0;e<mesh->Nelements;e++)
      invDegree[e] *= weight;

    //update with weight
    precon->o_invDegreeAP.copyFrom(invDegree);
  }
  free(invDegree);
}

void ellipticSetupSmootherDampedJacobi(elliptic_t *elliptic, precon_t *precon, 
                                       agmgLevel *level, dfloat lambda) {

  dfloat *invDiagA;
  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  ellipticBuildJacobi(elliptic,lambda, &invDiagA);

  precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
    
  level->device_smoother = dampedJacobi;

  //estimate the max eigenvalue of S*A
  dfloat rho = maxEigSmoothAx(elliptic, level);

  if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {

    level->smoother_params = (dfloat *) calloc(2,sizeof(dfloat));

    level->smoother_params[0] = rho;
    level->smoother_params[1] = rho/10.;

  } else {

    //set the stabilty weight (jacobi-type interation)
    dfloat weight = (4./3.)/rho;

    for (dlong n=0;n<mesh->Np*mesh->Nelements;n++)
      invDiagA[n] *= weight;

    //update diagonal with weight
    precon->o_invDiagA.copyFrom(invDiagA);
  }

  free(invDiagA);
}

static void eig(const int Nrows, double *A, double *WR, double *WI){

  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK  = (NB+2)*N;

  double *WORK  = new double[LWORK];
  double *VL  = new double[Nrows*Nrows];
  double *VR  = new double[Nrows*Nrows];

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
    VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);


  assert(INFO == 0);

  delete [] VL;
  delete [] VR;
  delete [] WORK;
}

dfloat maxEigSmoothAx(elliptic_t* elliptic, agmgLevel *level){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  const dlong N = level->Nrows;
  const dlong M = level->Ncols;

  int k = 10;

  hlong Nlocal = (hlong) level->Nrows;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat *Vx = (dfloat*) calloc(M, sizeof(dfloat));
  occa::memory *o_V = (occa::memory *) calloc(k+1, sizeof(occa::memory));
  
  occa::memory o_Vx  = mesh->device.malloc(M*sizeof(dfloat),Vx);
  occa::memory o_AVx = mesh->device.malloc(M*sizeof(dfloat),Vx);

  for(int i=0; i<=k; i++)
    o_V[i] = mesh->device.malloc(M*sizeof(dfloat),Vx);

  // generate a random vector for initial basis vector
  for (dlong i=0;i<N;i++) Vx[i] = (dfloat) drand48(); 

  //gather-scatter 
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    gsParallelGatherScatter(mesh->hostGsh, Vx, dfloatString, "add"); 
  
    for (dlong i=0;i<elliptic->Nmasked;i++) Vx[elliptic->maskIds[i]] = 0.;
  }

  o_Vx.copyFrom(Vx); //copy to device
  dfloat norm_vo = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_Vx, o_Vx);
  norm_vo = sqrt(norm_vo);

  ellipticScaledAdd(elliptic, 1./norm_vo, o_Vx, 0. , o_V[0]);

  for(int j=0; j<k; j++){
    // v[j+1] = invD*(A*v[j])
    level->device_Ax(level->AxArgs,o_V[j],o_AVx);
    level->device_smoother(level->smootherArgs, o_AVx, o_V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_V[i], o_V[j+1]);

      // v[j+1] = v[j+1] - hij*v[i]
      ellipticScaledAdd(elliptic, -hij, o_V[i], 1., o_V[j+1]);

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_V[j+1], o_V[j+1]);
      norm_vj = sqrt(norm_vj);
      ellipticScaledAdd(elliptic, 1/norm_vj, o_V[j+1], 0., o_V[j+1]);
      
      H[j+1+ j*k] = (double) norm_vj;
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  // free memory
  free(H);
  free(WR);
  free(WI);

  free(Vx);
  o_Vx.free();
  o_AVx.free();
  for(int i=0; i<=k; i++) o_V[i].free();
  free((void*)o_V);

  if((mesh->rank==0)&&(options.compareArgs("VERBOSE","TRUE"))) printf("weight = %g \n", rho);

  return rho;
}
