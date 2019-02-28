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

#include "adaptive.h"

/** Initialize the libraries that we are using
 *
 */
static int init_libs(MPI_Comm comm, int verbosity)
{
  int rank;
  ASD_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  int logpriority = ASD_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, logpriority);
  p4est_init(NULL, logpriority);

  int loglevel = ASD_MAX(ASD_LL_INFO - verbosity, ASD_LL_ALWAYS);
  asd_log_init(rank, stdout, loglevel);

  // add signal handler to get backtrace on abort
  asd_signal_handler_set();

  return loglevel;
}

static void print_precision()
{
  const char *comp = (sizeof(double) == sizeof(dfloat_t)) ? "double" : "single";
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("----- Precision ------------------------------------------");
  ASD_ROOT_INFO("compute precision = %s", comp);
  ASD_ROOT_INFO("----------------------------------------------------------");
}

int main(int argc, char **argv){
  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./adaptiveMain setupfile\n");

    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  int verbosity = 0;
  options.getArgs("VERBOSITY", verbosity);
  init_libs(comm, verbosity);
  print_precision();

  adaptive_t *adaptive = adaptiveSetup(options, comm);
  level_t *level = adaptive->lvl;

  dfloat lambda = 1;

  options.getArgs("LAMBDA", lambda);

  printf("lambda = %lf\n", lambda);
  
  dfloat *rhs   = (dfloat*) calloc(level->Np*level->Klocal, sizeof(dfloat));
  dfloat *x     = (dfloat*) calloc(level->Np*level->Klocal, sizeof(dfloat));
  dfloat *exact = (dfloat*) calloc(level->Np*level->Klocal, sizeof(dfloat));
  
  dfloat *vgeo = (dfloat*) calloc(NVGEO*level->Np*level->Klocal, sizeof(dfloat));
  dfloat *ggeo = (dfloat*) calloc(NGGEO*level->Np*level->Klocal, sizeof(dfloat));

  dfloat *vgeoGJ = (dfloat*) calloc(NVGEO*level->NqGJ*level->NqGJ*level->NqGJ*level->Klocal, sizeof(dfloat));
  dfloat *ggeoGJ = (dfloat*) calloc(NGGEO*level->NqGJ*level->NqGJ*level->NqGJ*level->Klocal, sizeof(dfloat));
  dfloat *gllw = (dfloat*) calloc(level->Nq, sizeof(dfloat));

  level->o_w.copyTo(gllw);

  level->o_vgeo.copyTo(vgeo, NVGEO*level->Np*level->Klocal*sizeof(dfloat), 0);
  level->o_ggeo.copyTo(ggeo, NGGEO*level->Np*level->Klocal*sizeof(dfloat), 0);

  level->o_vgeoGJ.copyTo(vgeoGJ, NVGEO*level->NqGJ*level->NqGJ*level->NqGJ*level->Klocal*sizeof(dfloat), 0);
  level->o_ggeoGJ.copyTo(ggeoGJ, NGGEO*level->NqGJ*level->NqGJ*level->NqGJ*level->Klocal*sizeof(dfloat), 0);

#if 0
  for(iint_t e=0;e<level->Klocal;++e){
    for(int k=0;k<level->Nq;++k){
      for(int j=0;j<level->Nq;++j){
	for(int i=0;i<level->Nq;++i){
	  int n = i + j*level->Nq + k*level->Nq*level->Nq;
	  dfloat J = vgeo[e*NVGEO*level->Np+n+level->Np*VGEO_J];
	  dfloat xn = vgeo[e*NVGEO*level->Np+n+level->Np*VGEO_X];
	  dfloat yn = vgeo[e*NVGEO*level->Np+n+level->Np*VGEO_Y];
	  dfloat zn = vgeo[e*NVGEO*level->Np+n+level->Np*VGEO_Z];
	  dfloat JWn = ggeo[e*NGGEO*level->Np+n+level->Np*GGEO_JW];
	  
	  dfloat mode = 2;
	  
	  iint_t id = n + e*level->Np;
	  rhs[id] =
	    JWn*(3*mode*mode*M_PI*M_PI+lambda)*cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);

	  exact[id] = cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);

	}
      }
    }
  }

  
#else

  
  int Nq = level->Nq;
  int Np = level->Np;
  int NqGJ = level->NqGJ;
  int NpGJ = NqGJ*NqGJ*NqGJ;
  dfloat bQQQ[NqGJ][NqGJ][NqGJ]; //  = (dfloat*) calloc(NpGJ, sizeof(dfloat));
  dfloat bQQN[NqGJ][NqGJ][NqGJ]; //  = (dfloat*) calloc(NpGJ, sizeof(dfloat));
  dfloat bQNN[NqGJ][NqGJ][NqGJ]; //  = (dfloat*) calloc(NpGJ, sizeof(dfloat));

  dfloat *IGJ = (dfloat*) calloc(NqGJ*Nq, sizeof(dfloat));

  level->o_IGJ.copyTo(IGJ);
  
  printf("IGJ=[\n");
  for(int i=0;i<NqGJ;++i){
    for(int a=0;a<Nq;++a){
      printf("%f ", IGJ[i*Nq+a]);
    }
    printf("\n");
  }
  printf("];\n");
  
  dfloat mode = 2;
  
  for(iint_t e=0;e<level->Klocal;++e){

    for(int k=0;k<NqGJ;++k){
      for(int j=0;j<NqGJ;++j){
	for(int i=0;i<NqGJ;++i){
	  int n = i + j*NqGJ + k*NqGJ*NqGJ;
	  dfloat JW = ggeoGJ[e*NGGEO*NpGJ+n+NpGJ*GGEO_JW];
	  dfloat xn = vgeoGJ[e*NVGEO*NpGJ+n+NpGJ*VGEO_X];
	  dfloat yn = vgeoGJ[e*NVGEO*NpGJ+n+NpGJ*VGEO_Y];
	  dfloat zn = vgeoGJ[e*NVGEO*NpGJ+n+NpGJ*VGEO_Z];

	  // dfloat JWn = ggeoGJ[e*NGGEO*NpGJ+n+NpGJ*GGEO_JW];
	  bQQQ[k][j][i] =
	    JW*(3*mode*mode*M_PI*M_PI+lambda)*cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);
	  //	  printf("JW = %lf, X=%lf,%lf,%lf bQQQ=%lf\n", JW, xn,yn,zn, bQQQ[k][j][i]);
	}
      }
    }
    
    for(int k=0;k<NqGJ;++k){
      for(int j=0;j<NqGJ;++j){
	for(int a=0;a<Nq;++a){
	  dfloat res = 0;
	  for(int i=0;i<NqGJ;++i){
	    res += IGJ[i*Nq+a]*bQQQ[k][j][i];
	  }
	  bQQN[k][j][a] = res;
	}
      }
    }

    for(int k=0;k<NqGJ;++k){
      for(int b=0;b<Nq;++b){
	for(int a=0;a<Nq;++a){
	  dfloat res = 0;
	  for(int j=0;j<NqGJ;++j){
	    res += IGJ[j*Nq+b]*bQQN[k][j][a];
	  }
	  bQNN[k][b][a] = res;
	}
      }
    }
  
    for(int c=0;c<Nq;++c){
      for(int b=0;b<Nq;++b){
	for(int a=0;a<Nq;++a){
	  dfloat res = 0;
	  for(int k=0;k<NqGJ;++k){
	    res += IGJ[k*Nq+c]*bQNN[k][b][a];
	  }

	  int id = a + b*Nq + c*Nq*Nq + e*level->Np;
	  rhs[id] = res;
  
	  int n = a + b*Nq + c*Nq*Nq;
	  
	  dfloat xn = vgeo[e*NVGEO*Np+n+Np*VGEO_X];
	  dfloat yn = vgeo[e*NVGEO*Np+n+Np*VGEO_Y];
	  dfloat zn = vgeo[e*NVGEO*Np+n+Np*VGEO_Z];
	  
	  exact[id] = cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);
	}
      }
    }
  }
#endif
  
  occa::memory o_rhs = adaptive->device.malloc(level->Np*level->Klocal*sizeof(dfloat), rhs);
  occa::memory o_x   = adaptive->device.malloc(level->Np*level->Klocal*sizeof(dfloat), x);

#if 0
  adaptive->dotMultiplyKernel(level->Np*level->Klocal, level->o_invDegree, o_rhs, o_rhs);

  adaptiveGatherScatter(level, o_rhs);

  o_rhs.copyTo(rhsb);
  
  for(int n=0;n<level->Np*level->Klocal;++n){
    if(fabs(rhs[n]-exact[n])>1e-12)
      printf("mismatch from %lf to %lf of %lf\n",
	     rhs[n], exact[n], fabs(rhs[n]-exact[n]));
  }
  MPI_Finalize();
  exit(0);
#endif
  
  dfloat tol = 1.e-6;
  adaptiveSolve(adaptive, lambda, tol, o_rhs, o_x);
  
  o_x.copyTo(x, level->Np*level->Klocal*sizeof(dfloat), 0);

  dfloat_t maxError = 0;
  for(iint_t n=0;n<level->Klocal*level->Np;++n){
    maxError = ASD_MAX(maxError, fabs(exact[n]-x[n]));
    x[n] -= exact[n];
  }
  printf("maxError = %le\n", maxError);
  
  o_x.copyFrom(x, level->Np*level->Klocal*sizeof(dfloat), 0);
  
  adaptivePlotVTUHex3D(adaptive, adaptive->lvl, 0, 0.0, "out", o_x);

  adaptive_free(adaptive);

  // close down MPI
  MPI_Finalize();

  return 0;
}
