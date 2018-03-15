#include <stdlib.h>
#include "mesh2D.h"

void meshGatherScatterQuad2D(mesh2D *mesh, dfloat *qL){

  for(int n=0;n<mesh->NgatherNodes;++n){
    const int Ngather = mesh->gatherCounts[n];
    const int offset  = mesh->gatherOffsets[n];

    dfloat sumqn = 0;
    for(int m=0;m<Ngather;++m){
      const int id = mesh->gatherIds[offset+m];
      sumqn += qL[id];
    }
    
    for(int m=0;m<Ngather;++m){
      const int id = mesh->gatherIds[offset+m];
      qL[id] = sumqn;
    }
  }
}
