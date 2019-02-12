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

typedef struct{

  dfloat VX;
  dfloat VY;

  dlong localId;
  hlong globalId;

}FEMverts_t;

typedef struct {

  dlong localId;
  hlong globalId;
  int ownerRank;

}parallelNode_t;

// compare on global owners
int parallelCompareOwnersAndGlobalId(const void *a, const void *b);

// compare on global indices
int parallelCompareGlobalId(const void *a, const void *b);

// compare xy coordinates
int parallelCompareFEMvertsLocation(const void *a, const void *b){
  dfloat NODETOL = 1e-6;

  FEMverts_t *fa = (FEMverts_t*) a;
  FEMverts_t *fb = (FEMverts_t*) b;

  if(fa->VX < fb->VX - NODETOL) return -1;
  if(fa->VX > fb->VX + NODETOL) return +1;

  if(fa->VY < fb->VY - NODETOL) return -1;
  if(fa->VY > fb->VY + NODETOL) return +1;

  return 0;
}

// compare local id
int parallelCompareFEMvertsLocalId(const void *a, const void *b){

  FEMverts_t *fa = (FEMverts_t*) a;
  FEMverts_t *fb = (FEMverts_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

int parallelCompareRowColumn(const void *a, const void *b);

void BuildFEMMatrixHex3D (mesh_t *femMesh, mesh_t *pmesh, dfloat lambda, dlong *localIds, hlong* globalNumbering,int *globalOwners,dlong *cnt, nonZero_t *A);



void adaptiveSEMFEMSetup(adaptive_t *adaptive, precon_t* precon, dfloat lambda) {

  setupAide options = adaptive->options;

  if (!(options.compareArgs("DISCRETIZATION", "CONTINUOUS"))) {
    printf("SEMFEM is supported for CONTINUOUS only\n");
    MPI_Barrier(adaptive->mesh->comm);
    MPI_Finalize();
    exit(0);
  }

  mesh_t* mesh = adaptive->mesh; //original mesh

  //  mesh_t* pmesh = (mesh_t*) calloc (1,sizeof(mesh_t)); //partially assembled fem mesh (result of projecting sem element to larger space)
  mesh_t* pmesh = new mesh_t[1];
  
  //  precon->femMesh = (mesh_t*) calloc (1,sizeof(mesh_t)); //full fem mesh
  precon->femMesh = new mesh_t[1];
  
  mesh_t *femMesh = precon->femMesh;

  memcpy(pmesh  ,mesh,sizeof(mesh_t));
  memcpy(femMesh,mesh,sizeof(mesh_t));


  //now build the full degree 1 fem mesh
  int femN = 1; //degree of fem approximation

  /* allocate space for node coordinates */
  femMesh->Nelements = mesh->NelFEM*mesh->Nelements;
  femMesh->EToV = (hlong*) calloc(femMesh->Nelements*femMesh->Nverts, sizeof(hlong));
  femMesh->EX = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  femMesh->EY = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));
  if (adaptive->dim==3)
    femMesh->EZ = (dfloat*) calloc(femMesh->Nverts*femMesh->Nelements, sizeof(dfloat));

  dlong *localIds = (dlong *) calloc(femMesh->Nverts*femMesh->Nelements,sizeof(dlong));

  // dlong NFEMverts = mesh->Nelements*mesh->NpFEM;
  for(dlong e=0;e<mesh->Nelements;++e){
    for (int n=0;n<mesh->NelFEM;n++) {
      dlong id[femMesh->Nverts];

      dlong femId = e*mesh->NelFEM*mesh->Nverts+n*mesh->Nverts;

      for (int i=0;i<femMesh->Nverts;i++) {
        //local ids in the subelement fem grid
        id[i] = e*mesh->NpFEM + mesh->FEMEToV[n*mesh->Nverts+i];

        /* read vertex triplet for triangle */
        femMesh->EToV[femId+i] = pmesh->globalIds[id[i]];

        femMesh->EX[femId+i] = pmesh->x[id[i]];
        femMesh->EY[femId+i] = pmesh->y[id[i]];
        if (adaptive->dim==3)
          femMesh->EZ[femId+i] = pmesh->z[id[i]];

      }

      localIds[femId+0] = id[0];
      localIds[femId+1] = id[1];
      localIds[femId+2] = id[3];  //need to swap this as the Np nodes are ordered [0,1,3,2,4,5,7,6] in a degree 1 element
      localIds[femId+3] = id[2];
      localIds[femId+4] = id[4];
      localIds[femId+5] = id[5];
      localIds[femId+6] = id[7];
      localIds[femId+7] = id[6];
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(femMesh);

  meshLoadReferenceNodesHex3D(femMesh, femN);

  int *faceFlag = (int*) calloc(pmesh->Np*pmesh->Nfaces,sizeof(int));
  for (int f=0;f<pmesh->Nfaces;f++) {
    for (int n=0;n<pmesh->Nfp;n++) {
      int id = pmesh->faceNodes[f*pmesh->Nfp+n];
      faceFlag[f*pmesh->Np + id] = 1; //flag the nodes on this face
    }
  }

  //map from faces of fem sub-elements to the macro element face number
  int *femFaceMap = (int*) calloc(mesh->NelFEM*femMesh->Nfaces,sizeof(int));
  for (int n=0;n<mesh->NelFEM*femMesh->Nfaces;n++) femFaceMap[n] = -1;

  for (int n=0;n<mesh->NelFEM;n++) {
    for (int f=0;f<femMesh->Nfaces;f++) {

      for (int face=0; face<pmesh->Nfaces;face++) {

        //count the nodes on this face which are on a macro face
        int NvertsOnFace = 0;
        for (int i=0;i<femMesh->Nfp;i++){
          int id = femMesh->faceNodes[f*femMesh->Nfp+i];
          int v  = mesh->FEMEToV[n*pmesh->Nverts+id];
          NvertsOnFace += faceFlag[face*pmesh->Np + v];
        }
        if (NvertsOnFace == femMesh->Nfp)
          femFaceMap[n*femMesh->Nfaces+f] = face; //on macro face
      }
    }
  }

  //fill the boundary flag array
  femMesh->EToB = (int*) calloc(femMesh->Nelements*femMesh->Nfaces, sizeof(int));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->NelFEM;n++) {
      for (int f=0;f<femMesh->Nfaces;f++) {
        int face = femFaceMap[n*femMesh->Nfaces+f];
        if (face>-1) {
          femMesh->EToB[(e*mesh->NelFEM +n)*femMesh->Nfaces +f] = mesh->EToB[e*mesh->Nfaces + face];
        }
      }
    }
  }
  free(faceFlag);
  free(femFaceMap);

  meshPhysicalNodesHex3D(femMesh);
  meshGeometricFactorsHex3D(femMesh);
  meshHaloSetup(femMesh);
  meshConnectFaceNodes3D(femMesh);
  meshSurfaceGeometricFactorsHex3D(femMesh);

  // global nodes
  meshParallelConnectNodes(femMesh);

  dlong Ntotal = pmesh->Np*pmesh->Nelements;
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;

  pmesh->maskedGlobalIds = (hlong *) calloc(Ntotal,sizeof(hlong));
  memcpy(pmesh->maskedGlobalIds, pmesh->globalIds, Ntotal*sizeof(hlong));

  //mask using the original mask
  for (dlong n=0;n<adaptive->Nmasked;n++)
    pmesh->maskedGlobalIds[adaptive->maskIds[n]] = 0;

  //build masked gs handle
  precon->FEMogs = ogsSetup(Ntotal, pmesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = precon->FEMogs->Ngather;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc(Ngather,sizeof(hlong));
  int   *owner     = (int *) calloc(Ngather,sizeof(int));

  // every gathered degree of freedom has its own global id
  hlong *globalStarts = (hlong *) calloc(mesh->size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<precon->FEMogs->Ngather;n++) {
    globalIds[n] = n + globalStarts[mesh->rank];
    owner[n] = mesh->rank;
  }

  //scatter this numbering to the original nodes
  hlong *globalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  int *globalOwners = (int *) calloc(Ntotal,sizeof(int));
  for (dlong n=0;n<Ntotal;n++) globalNumbering[n] = -1;
  ogsScatter(globalNumbering, globalIds, ogsHlong, ogsAdd, precon->FEMogs);
  ogsScatter(globalOwners, owner, ogsInt, ogsAdd, precon->FEMogs);


  free(globalIds); free(owner);

  if (mesh->rank==0) printf("Building full SEMFEM matrix..."); fflush(stdout);

  // Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = femMesh->Np*femMesh->Np*femMesh->Nelements;

  dlong cnt =0;
  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(mesh->size, sizeof(int));
  int *AsendOffsets = (int*) calloc(mesh->size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(mesh->size+1, sizeof(int));

  //Build unassembed non-zeros
  BuildFEMMatrixHex3D(femMesh,pmesh,lambda, localIds, globalNumbering, globalOwners,&cnt,sendNonZeros); 

  // Make the MPI_NONZERO_T data type
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[4] = {MPI_HLONG, MPI_HLONG, MPI_INT, MPI_DFLOAT};
  int blength[4] = {1, 1, 1, 1};
  MPI_Aint addr[4], displ[4];
  MPI_Get_address ( &(sendNonZeros[0]          ), addr+0);
  MPI_Get_address ( &(sendNonZeros[0].col      ), addr+1);
  MPI_Get_address ( &(sendNonZeros[0].ownerRank), addr+2);
  MPI_Get_address ( &(sendNonZeros[0].val      ), addr+3);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  MPI_Type_create_struct (4, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  // count how many non-zeros to send to each process
  for(dlong n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank]++;

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh->comm);

  // find send and recv offsets for gather
  dlong nnz = 0;
  for(int r=0;r<mesh->size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    nnz += ArecvCounts[r];
  }

  nonZero_t *A = (nonZero_t*) calloc(nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_NONZERO_T,
		A, ArecvCounts, ArecvOffsets, MPI_NONZERO_T,
		mesh->comm);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort(A, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<nnz;++n){
    if(A[n].row == A[cnt].row && A[n].col == A[cnt].col){
      A[cnt].val += A[n].val;
    } else{
      ++cnt;
      A[cnt] = A[n];
    }
  }
  if (nnz) cnt++;
  nnz = cnt;

  if(mesh->rank==0) printf("done.\n");

  MPI_Barrier(mesh->comm);
  MPI_Type_free(&MPI_NONZERO_T);

  hlong *Rows = (hlong *) calloc(nnz, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnz, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

  for (dlong n=0;n<nnz;n++) {
    Rows[n] = A[n].row;
    Cols[n] = A[n].col;
    Vals[n] = A[n].val;
  }
  free(A);

  precon->parAlmond = parAlmond::Init(mesh->device, mesh->comm, options);
  parAlmond::AMGSetup(precon->parAlmond,
                     globalStarts,
                     nnz,
                     Rows,
                     Cols,
                     Vals,
                     adaptive->allNeumann,
                     adaptive->allNeumannPenalty);
  free(Rows); free(Cols); free(Vals);

  if (options.compareArgs("VERBOSE", "TRUE"))
      parAlmond::Report(precon->parAlmond);

  parAlmond::multigridLevel *baseLevel = precon->parAlmond->levels[0];
  precon->rhsG = (dfloat*) calloc(baseLevel->Ncols,sizeof(dfloat));
  precon->xG   = (dfloat*) calloc(baseLevel->Ncols,sizeof(dfloat));
  precon->o_rhsG = mesh->device.malloc(baseLevel->Ncols*sizeof(dfloat));
  precon->o_xG   = mesh->device.malloc(baseLevel->Ncols*sizeof(dfloat));
}

void BuildFEMMatrixHex3D(mesh_t *femMesh, mesh_t *pmesh, dfloat lambda,
                        dlong *localIds, hlong* globalNumbering, int *globalOwners,
                        dlong *cnt, nonZero_t *A) {

#pragma omp parallel for
  for (dlong e=0;e<femMesh->Nelements;e++) {
    for (int nz=0;nz<femMesh->Nq;nz++) {
      for (int ny=0;ny<femMesh->Nq;ny++) {
	for (int nx=0;nx<femMesh->Nq;nx++) {
	  dlong nn = nx+ny*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
	  dlong idn = localIds[e*femMesh->Np + nn];
	  if (globalNumbering[idn]<0) continue; //skip masked nodes

	  for (int mz=0;mz<femMesh->Nq;mz++) {
	    for (int my=0;my<femMesh->Nq;my++) {
	      for (int mx=0;mx<femMesh->Nq;mx++) {
		dlong mm = mx+my*femMesh->Nq+mz*femMesh->Nq*femMesh->Nq;
		dlong idm = localIds[e*femMesh->Np + mm];
		if (globalNumbering[idm]<0) continue; //skip masked nodes

		int id;
		dfloat val = 0.;

		if ((ny==my)&&(nz==mz)) {
		  for (int k=0;k<femMesh->Nq;k++) {
		    id = k+ny*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		    dfloat Grr = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G00ID*femMesh->Np];

		    val += Grr*femMesh->D[nx+k*femMesh->Nq]*femMesh->D[mx+k*femMesh->Nq];
		  }
		}

		if (nz==mz) {
		  id = mx+ny*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		  dfloat Grs = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G01ID*femMesh->Np];
		  val += Grs*femMesh->D[nx+mx*femMesh->Nq]*femMesh->D[my+ny*femMesh->Nq];

		  id = nx+my*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		  dfloat Gsr = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G01ID*femMesh->Np];
		  val += Gsr*femMesh->D[mx+nx*femMesh->Nq]*femMesh->D[ny+my*femMesh->Nq];
		}

		if (ny==my) {
		  id = mx+ny*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		  dfloat Grt = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G02ID*femMesh->Np];
		  val += Grt*femMesh->D[nx+mx*femMesh->Nq]*femMesh->D[mz+nz*femMesh->Nq];

		  id = nx+ny*femMesh->Nq+mz*femMesh->Nq*femMesh->Nq;
		  dfloat Gst = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G02ID*femMesh->Np];
		  val += Gst*femMesh->D[mx+nx*femMesh->Nq]*femMesh->D[nz+mz*femMesh->Nq];
		}

		if ((nx==mx)&&(nz==mz)) {
		  for (int k=0;k<femMesh->Nq;k++) {
		    id = nx+k*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		    dfloat Gss = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G11ID*femMesh->Np];

		    val += Gss*femMesh->D[ny+k*femMesh->Nq]*femMesh->D[my+k*femMesh->Nq];
		  }
		}

		if (nx==mx) {
		  id = nx+my*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		  dfloat Gst = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G12ID*femMesh->Np];
		  val += Gst*femMesh->D[ny+my*femMesh->Nq]*femMesh->D[mz+nz*femMesh->Nq];

		  id = nx+ny*femMesh->Nq+mz*femMesh->Nq*femMesh->Nq;
		  dfloat Gts = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G12ID*femMesh->Np];
		  val += Gts*femMesh->D[my+ny*femMesh->Nq]*femMesh->D[nz+mz*femMesh->Nq];
		}

		if ((nx==mx)&&(ny==my)) {
		  for (int k=0;k<femMesh->Nq;k++) {
		    id = nx+ny*femMesh->Nq+k*femMesh->Nq*femMesh->Nq;
		    dfloat Gtt = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + G22ID*femMesh->Np];

		    val += Gtt*femMesh->D[nz+k*femMesh->Nq]*femMesh->D[mz+k*femMesh->Nq];
		  }
		}

		if ((nx==mx)&&(ny==my)&&(nz==mz)) {
		  id = nx + ny*femMesh->Nq+nz*femMesh->Nq*femMesh->Nq;
		  dfloat JW = femMesh->ggeo[e*femMesh->Np*femMesh->Nggeo + id + GWJID*femMesh->Np];
		  val += JW*lambda;
		}

		// pack non-zero
		dfloat nonZeroThreshold = 1e-7;
		if (fabs(val) >= nonZeroThreshold) {
#pragma omp critical
		  {
		    A[*cnt].val = val;
		    A[*cnt].row = globalNumbering[idn];
		    A[*cnt].col = globalNumbering[idm];
		    A[*cnt].ownerRank = globalOwners[idn];
		    (*cnt)++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
