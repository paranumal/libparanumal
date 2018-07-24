/* compile with C compiler (not C++) */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"

#include "types.h"

void *gsParallelGatherScatterSetup(MPI_Comm meshComm,
				   dlong NuniqueBases,
				   hlong *gatherGlobalNodes,
                                   int verbose){

  /* gslib stuff */
  comm_ext world;
  struct comm com;
  
  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)meshComm; // MPI_COMM_WORLD;
  
  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);
  
  dlong n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = (slong) gatherGlobalNodes[n];
  }

  struct gs_data *gsh = gs_setup(id, NuniqueBases, &com, 0, gs_auto, verbose);

  free(id);

  return gsh;
}


void gsParallelGatherScatterDestroy(void *gsh){

  gs_free(gsh);

}
