/* compile with C compiler (not C++) */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#define UNDERSCORE 1
#define USE_NAIVE_BLAS 
#define NO_NEX_EXITT 1
#define GLOBAL_LONG_LONG 1
#define PREFIX jl_

#define MPI 1

#include "gslib.h"

void *meshParallelGatherScatterSetup(int NuniqueBases,
				     int *gatherGlobalNodes){
  
  /* gslib stuff */
  comm_ext world;
  struct comm com;
  
  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)MPI_COMM_WORLD;
  
  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);
  
  int n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = gatherGlobalNodes[n];
  }

  struct gs_data *gsh = gs_setup(id, NuniqueBases, &com, 0, gs_auto, 1);

  free(id);

  return gsh;
}
