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

void serialTic();
double serialElapsed();

// Scalable Non-blocking Preconditioned Conjugate Gradient Methods, by Paul R. Eller && William Gropp
// http://eller3.web.engr.illinois.edu/eller_sc16_presentation.pdf

int nbpcg(elliptic_t* elliptic, dfloat lambda, 
	  occa::memory &o_r, occa::memory &o_x, 
	  const dfloat tol, const int MAXIT){
  
  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;
  
  const cgOptions_t &cgOptions = elliptic->cgOptions;

  int fixedIterationCountFlag = 0;
  if(options.compareArgs("FIXED ITERATION COUNT", "TRUE")){
    fixedIterationCountFlag = 1;
  }
  
  // register scalars
  dfloat normz0 = 0;
  dfloat zdotz0 = 0;

  dfloat alpha0 = 0;
  dfloat beta0  = 0;
  dfloat gamma0 = 0;
  dfloat delta0 = 0;

  dfloat gamma1 = 0; // history gamma
  
  dfloat one = 1;
  
  double serialElapsedReduction = 0, serialElapsedAx = 0, serialElapsedGatherScatter = 0;

  double pcgStart = MPI_Wtime();
  
  /*aux variables */
  occa::memory &o_p  = elliptic->o_pcgWork[0];
  occa::memory &o_s  = elliptic->o_pcgWork[1];
  occa::memory &o_S  = elliptic->o_pcgWork[2];
  occa::memory &o_z  = elliptic->o_pcgWork[3];
  occa::memory &o_Z  = elliptic->o_pcgWork[4];
  occa::memory &o_Ax  = elliptic->o_pcgWork[5];

  MPI_Request request;
  MPI_Status  status;

  dfloat *localdots = (dfloat*) calloc(2, sizeof(dfloat));
  dfloat *globaldots = (dfloat*) calloc(2, sizeof(dfloat));
  
  // Ax = A*x
  ellipticOperator(elliptic, lambda, o_x, o_Ax, dfloatString); // WRONG FOR IPDG
  
  // subtract r = b - A*x
  ellipticScaledAdd(elliptic, -one, o_Ax, one, o_r);

  // z = M*r [ Gropp notation ]
  ellipticPreconditioner(elliptic, lambda, o_r, o_z);
  
  // set alpha = 0 to get
  // r.z and z.z
  // [
  alpha0 = 0;
  ellipticNonBlockingUpdate2NBPCG(elliptic, o_s, o_S, alpha0, o_r, o_z, localdots, globaldots, &request);

  ellipticOperator(elliptic, lambda, o_z, o_Z, dfloatString);

  MPI_Wait(&request, &status);
  gamma0 = globaldots[0]; // rdotz
  zdotz0 = globaldots[1];
  normz0 = sqrt(zdotz0);
  // ]

  dfloat TOL = mymax(gamma0*tol*tol, tol*tol);

  int iter;

  beta0 = 0;
  for(iter=1;iter<=MAXIT;++iter){

    // p <= z + beta*p
    // s <= Z + beta*s
    // delta <= pdots
    ellipticNonBlockingUpdate1NBPCG(elliptic, o_z, o_Z, beta0, o_p, o_s, localdots, globaldots, &request);

    // S = M*s
    ellipticPreconditioner(elliptic, lambda, o_s, o_S);

    // block for delta
    MPI_Wait(&request, &status);
    delta0 = globaldots[0];

    // alpha = gamma/delta
    alpha0 = gamma0/delta0;

    // r <= r - alpha*s    
    // z <= z - alpha*S
    // r.z
    // z.z
    ellipticNonBlockingUpdate2NBPCG(elliptic, o_s, o_S, alpha0, o_r, o_z, localdots, globaldots, &request);

    // x <= x + alpha*p (delayed)
    ellipticScaledAdd(elliptic, alpha0, o_p, one, o_x);
    
    // Z = A*z
    ellipticOperator(elliptic, lambda, o_z, o_Z, dfloatString); // WRONG FOR IPDG        
    
    // block for delta
    MPI_Wait(&request, &status);
    gamma1 = gamma0;
    gamma0 = globaldots[0]; // gamma = r.z
    zdotz0 = globaldots[1]; // 
    normz0 = sqrt(zdotz0);

    beta0 = gamma0/gamma1;
    
    if (cgOptions.verbose&&(mesh->rank==0)) {

      if(gamma0<0)
	printf("WARNING CG: zdotz = %17.15lf\n", zdotz0);
      
      printf("CG: it %d z norm %12.12le gamma = %le zdotz = %le\n",
	     iter, normz0, gamma0, zdotz0);
    }

    if(gamma0<=TOL && !fixedIterationCountFlag) break;
    
  }

  free(localdots);
  free(globaldots);
  
  return iter;
}


