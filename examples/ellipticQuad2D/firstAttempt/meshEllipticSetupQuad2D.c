#include "mpi.h"
#include "mesh2D.h"

void meshEllipticSetupQuad2D(mesh2D *mesh){

  // set up global numbering of nodes
  meshNumberNodes2D(mesh);
  
  // set up gather/scatter ops (serial)
  void meshParallelScatterSetupQuad2D(mesh2D *mesh);
  void meshParallelGatherSetupQuad2D(mesh2D *mesh);
  
  meshParallelGatherSetupQuad2D(mesh); 
  meshParallelScatterSetupQuad2D(mesh);

  // need to set up mask and multiplicity ?
}
