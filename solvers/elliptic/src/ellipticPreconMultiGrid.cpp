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

#include "ellipticPrecon.hpp"

void interpolateFieldQuad2D(mesh_t &mesh, int outNq, dfloat *outr, dfloat *outs, dfloat *inx, dfloat *outx){
  int inNq = mesh.Nq;
  dfloat tmpxOI[outNq][inNq];

  dfloat I[outNq][inNq];
  mesh.InterpolationMatrix1D(mesh.N,
			     inNq, mesh.r,
			     outNq, outr,
			     I[0]);
  
  for(dlong e=0;e<mesh.Nelements;++e){
    // interpolate in 's'
    for(int i=0;i<inNq;++i){
      for(int j=0;j<outNq;++j){
	tmpxOI[j][i] = 0;
	for(int n=0;n<inNq;++n){
	  tmpxOI[j][i] += I[j][n]*inx[e*mesh.Np + n*inNq + i];
	}
      }
    }
    
    // interpolate in 'r'
    for(int j=0;j<outNq;++j){
      for(int i=0;i<outNq;++i){
	dfloat xOO = 0;
	for(int n=0;n<inNq;++n){
	  xOO += I[i][n]*tmpxOI[j][n];
	}
	outx[e*outNq*outNq + j*outNq + i] = xOO;
      }
    }
  }
}


void interpolateFieldHex3D(mesh_t &mesh, int outNq, dfloat *outr, dfloat *outs, dfloat *outt, dfloat *inx, dfloat *outx){
  int inNq = mesh.Nq;
  dfloat tmpxOII[outNq][inNq][inNq];
  dfloat tmpxOOI[outNq][outNq][inNq];

  dfloat I[outNq][inNq];
  mesh.InterpolationMatrix1D(mesh.N,
			     inNq, mesh.r,
			     outNq, outr,
			     I[0]);

  for(dlong e=0;e<mesh.Nelements;++e){
    // interpolate in 't'
    for(int j=0;j<inNq;++j){
      for(int i=0;i<inNq;++i){
	for(int k=0;k<outNq;++k){
	  tmpxOII[k][j][i] = 0;
	  for(int n=0;n<inNq;++n){
	    tmpxOII[k][j][i] += I[k][n]*inx[e*mesh.Np + n*inNq*inNq + j*inNq + i];
	  }
	}
      }
    }
    // interpolate in 's'
    for(int k=0;k<outNq;++k){
      for(int i=0;i<inNq;++i){
	for(int j=0;j<outNq;++j){
	  tmpxOOI[k][j][i] = 0;
	  for(int n=0;n<inNq;++n){
	    tmpxOOI[k][j][i] += I[j][n]*tmpxOII[k][n][i];
	  }
	}
      }
    }
    
    // interpolate in 'r'
    for(int k=0;k<outNq;++k){
      for(int j=0;j<outNq;++j){
	for(int i=0;i<outNq;++i){
	  dfloat xOOO = 0;
	  for(int n=0;n<inNq;++n){
	    xOOO += I[i][n]*tmpxOOI[k][j][n];
	  }
	  outx[e*outNq*outNq*outNq + k*outNq*outNq + j*outNq + i] = xOOO;
	}
      }
    }
  }
}
    

void interpolateField(mesh_t &mesh, int outNq, dfloat *outr, dfloat *outs, dfloat *outt, dfloat *inx, dfloat *outx){

  if(mesh.elementType==HEXAHEDRA)
    interpolateFieldHex3D(mesh, outNq, outr, outs, outt, inx, outx);

  if(mesh.elementType==QUADRILATERALS)
    interpolateFieldQuad2D(mesh, outNq, outr, outs, inx, outx);
  
}

void interpolatePhysicalNodes(mesh_t &meshI, mesh_t &meshO){

  printf("interpolating physical nodes from N=%d to %d\n", meshI.N, meshO.N);
  
  // interpolate coordinates from degree N mesh to new mesh
  if(1){
    interpolateField(meshI, meshO.Nq, meshO.r, meshO.s, meshO.t, meshI.x, meshO.x);
    interpolateField(meshI, meshO.Nq, meshO.r, meshO.s, meshO.t, meshI.y, meshO.y);
    if(meshI.dim==3)
      interpolateField(meshI, meshO.Nq, meshO.r, meshO.s, meshO.t, meshI.z, meshO.z);
  }
  // compute geometric factors
  meshO.GeometricFactors();
  
  // compute surface geofacs
  meshO.SurfaceGeometricFactors();
  
  // reup geofacs
  meshO.o_ggeo.copyFrom(meshO.ggeo);
  meshO.o_sgeo.copyFrom(meshO.sgeo);
  meshO.o_vgeo.copyFrom(meshO.vgeo);

  // initialize cubature
  //meshO.CubatureSetup(meshI.cubN, meshI.cubatureType);
  meshO.CubatureSetup(meshO.N+1, "GL"); // meshI.cubatureType); // was N
}


// Matrix-free p-Multigrid levels followed by AMG
void MultiGridPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  //just pass to parAlmond
  parAlmond.Operator(o_r, o_Mr);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MultiGridPrecon::MultiGridPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, mesh.comm) {

  int Nf = mesh.N;
  int Nc = Nf;
  int NpFine   = mesh.Np;
  int NpCoarse = mesh.Np;

  MGLevel* prevLevel=nullptr;
  MGLevel* currLevel=nullptr;

  while(Nc>1) {
    //build mesh and elliptic objects for this degree
    mesh_t &meshF = mesh.SetupNewDegree(Nf);

    if(mesh.elementType==HEXAHEDRA || mesh.elementType==QUADRILATERALS){
      
      if(meshF.N!=mesh.N){
	interpolatePhysicalNodes(mesh, meshF);
      }
    }

    elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

    //    void plotGeometry(elliptic_t &elliptic, int num);
    //    plotGeometry(ellipticF, Nc);
    
    //share masking data with previous MG level
    if (prevLevel) {
      prevLevel->meshC = &meshF;
      prevLevel->ogsMaskedC = ellipticF.ogsMasked;
    }

    //find the degree of the next level
    if (settings.compareSetting("MULTIGRID COARSENING","ALLDEGREES")) {
      Nc = Nf-1;
    } else if (settings.compareSetting("MULTIGRID COARSENING","HALFDEGREES")) {
      Nc = mymax(1,(Nf+1)/2);
    } else { //default "HALFDOFS"
      // pick the degrees so the dofs of each level halfs (roughly)
      while (NpCoarse > NpFine/2 && Nc>1) {
        Nc--;
        switch(mesh.elementType){
          case TRIANGLES:
            NpCoarse = ((Nc+1)*(Nc+2))/2; break;
          case QUADRILATERALS:
            NpCoarse = (Nc+1)*(Nc+1); break;
          case TETRAHEDRA:
            NpCoarse = ((Nc+1)*(Nc+2)*(Nc+3))/6; break;
          case HEXAHEDRA:
            NpCoarse = (Nc+1)*(Nc+1)*(Nc+1); break;
        }
      }
    }

    //set Npcoarse
    switch(mesh.elementType){
      case TRIANGLES:
        NpCoarse = ((Nc+1)*(Nc+2))/2; break;
      case QUADRILATERALS:
        NpCoarse = (Nc+1)*(Nc+1); break;
      case TETRAHEDRA:
        NpCoarse = ((Nc+1)*(Nc+2)*(Nc+3))/6; break;
      case HEXAHEDRA:
        NpCoarse = (Nc+1)*(Nc+1)*(Nc+1); break;
    }

    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = ellipticF.ogsMasked->Ngather;
      Ncols = Nrows + ellipticF.ogsMasked->NgatherHalo;
    } else {
      Nrows = meshF.Nelements*meshF.Np;
      Ncols = Nrows + meshF.totalHaloPairs*mesh.Np;
    }

    //make a multigrid level
    currLevel = new MGLevel(ellipticF, Nrows, Ncols, Nc, NpCoarse);
    parAlmond.AddLevel(currLevel);

    Nf = Nc;
    NpFine = NpCoarse;
    prevLevel = currLevel;
  }

  //build matrix at degree 1
  mesh_t &meshF = mesh.SetupNewDegree(1);

 
  // reset geometry and quadrature to match finest grid
  if(mesh.elementType==HEXAHEDRA || mesh.elementType==QUADRILATERALS){

    if(meshF.N!=mesh.N){
      //    interpolatePhysicalNodes(mesh, meshF);
      meshF.CubatureSetup(meshF.N, "GLL"); // was N
    }else{
      meshF.CubatureSetup(mesh.N, "GLL"); // was N
    }
  }

  
  elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

  void plotGeometry(elliptic_t &elliptic, int num);
  plotGeometry(ellipticF, meshF.N);
  
  //share masking data with previous MG level
  if (prevLevel) {
    prevLevel->meshC = &meshF;
    prevLevel->ogsMaskedC = ellipticF.ogsMasked;
  }

  //build full A matrix and pass to parAlmond
  parAlmond::parCOO A(elliptic.platform, mesh.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticF.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticF.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = mesh.rank;
  int size = mesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = (dlong) (A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);
  free(null);

  //report
  parAlmond.Report();
}
