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

  adaptive_t *adaptive = adaptive_new(options, comm);
  level_t *level = adaptive->lvl;

  dfloat *b = (dfloat*) calloc(level->Np*level->Klocal, sizeof(dfloat));
  dfloat *x = (dfloat*) calloc(level->Np*level->Klocal, sizeof(dfloat));

  for(iint_t e=0;e<level->Klocal;++e){
    for(int n=0;n<level->Np;++n){
      b[e*level->Np+n] = drand48();
    }
  }
  
  occa::memory o_b = adaptive->device.malloc(level->Np*level->Klocal*sizeof(dfloat), b);
  occa::memory o_x = adaptive->device.malloc(level->Np*level->Klocal*sizeof(dfloat), x);

  dfloat lambda = 10, tol = 1.e-6;
  adaptiveSolve(adaptive, lambda, tol, o_b, o_x);
  
  adaptive_free(adaptive);

  // close down MPI
  MPI_Finalize();

  return 0;
}
