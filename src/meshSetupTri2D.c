#include "mesh2D.h"

mesh2D *meshSetupTri2D(char *filename, int N){

  // read chunk of elements
  mesh2D *mesh = meshParallelReaderTri2D(filename);

  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition2D(mesh);

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // load reference (r,s) element nodes
  meshLoadReferenceNodesTri2D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTri2D(mesh);

  // compute geometric factors
  meshGeometricFactorsTri2D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTri2D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  // initialize LSERK4 time stepping coefficients
  int Nrk = 5;

  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  // added one more for advanced time step
  dfloat rkc[6] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0 ,
		   1.0};

  mesh->Nrk = Nrk;
  memcpy(mesh->rka, rka, Nrk*sizeof(dfloat));
  memcpy(mesh->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(mesh->rkc, rkc, (Nrk+1)*sizeof(dfloat));

  //Adam-Bashforth
  mesh->mrab[0] = 23./12.;
  mesh->mrab[1] = -4./3.;
  mesh->mrab[2] =  5./12.;

  //AB half step
  mesh->mrabb[0] = 17./24.;
  mesh->mrabb[1] = -7./24.;
  mesh->mrabb[2] =  2./24.;

  int Nimex = 4;


  dfloat ImB[4] ={0.0,
              673488652607.0 /2334033219546.0,
              493801219040.0/853653026979.0,
              184814777513.0/1389668723319.0 };


  dfloat ImC[4] = { 0.0,
                    3375509829940.0/4525919076317.0,
                    272778623835.0/1039454778728.0,
                    1.0};



  dfloat ImAd[4] = {0.0,
                3375509829940.0/4525919076317.0,
                566138307881.0/912153721139.0,
                184814777513.0/1389668723319.0};

  dfloat ImAmBim[4] = {0.0,
                    0.0,
                  -11712383888607531889907.0/32694570495602105556248.0 - 673488652607.0 /2334033219546.0,
                   0.0};

  dfloat ImAmBex[4] = {0.0,
                    3375509829940.0/4525919076317.0,
                    272778623835.0/1039454778728.0 - 673488652607.0 /2334033219546.0,
                    1660544566939.0/2334033219546.0-493801219040.0/853653026979.0 };


  mesh->Nimex = Nimex;
  memcpy(mesh->LsimexB, ImB, Nimex*sizeof(dfloat));
  memcpy(mesh->LsimexC, ImC, Nimex*sizeof(dfloat));
  memcpy(mesh->LsimexAd, ImAd, Nimex*sizeof(dfloat));
  memcpy(mesh->LsimexABi, ImAmBim, Nimex*sizeof(dfloat));
  memcpy(mesh->LsimexABe, ImAmBex, Nimex*sizeof(dfloat));

  // Clasical Adams-Bashforth Coefficients


  return mesh;
}
