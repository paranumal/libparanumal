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

#include "elliptic.hpp"
#include "mesh/meshDefines2D.h"
#include "mesh/meshDefines3D.h"

// compare on global indices
int parallelCompareRowColumn(const void *a, const void *b){

  parAlmond::parCOO::nonZero_t *fa = (parAlmond::parCOO::nonZero_t*) a;
  parAlmond::parCOO::nonZero_t *fb = (parAlmond::parCOO::nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;
}

void compressMatrix(mesh_t &mesh, ogs_t *ogsMasked,
		    parAlmond::parCOO::nonZero_t *AL,
		    dlong cnt,
		    parAlmond::parCOO& A){

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked->Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  A.globalColStarts = (hlong*) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, A.globalRowStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }
  
  int *AsendCounts  = (int*) calloc(mesh.size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(mesh.size, sizeof(int));
  int *AsendOffsets = (int*) calloc(mesh.size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(mesh.size+1, sizeof(int));
  
  // sort by row ordering
  qsort(AL, cnt, sizeof(parAlmond::parCOO::nonZero_t), parallelCompareRowColumn);

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = AL[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, mesh.comm);

  // find send and recv offsets for gather
  A.nnz = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries = (parAlmond::parCOO::nonZero_t*) calloc(A.nnz, sizeof(parAlmond::parCOO::nonZero_t));

  // determine number to receive
  MPI_Alltoallv(AL, AsendCounts, AsendOffsets, parAlmond::MPI_NONZERO_T,
		A.entries, ArecvCounts, ArecvOffsets, parAlmond::MPI_NONZERO_T,
		mesh.comm);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((A.entries), A.nnz, sizeof(parAlmond::parCOO::nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
      A.entries[cnt].val += A.entries[n].val;
    }
    else{
      double tol = 1e-12;
      if(fabs(A.entries[n].val)>tol){
	++cnt;
	A.entries[cnt] = A.entries[n];
      }
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

  if(mesh.rank==0) printf("done %d entries.\n", cnt);

  MPI_Barrier(mesh.comm);
  free(AL);
  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

}

void elliptic_t::BuildOperatorMatrixContinuous(parAlmond::parCOO& A) {

  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;
  parAlmond::parCOO::nonZero_t *AL =
    (parAlmond::parCOO::nonZero_t*) calloc(nnzLocal, sizeof(parAlmond::parCOO::nonZero_t));

  if(mesh.rank==0) {printf("Building full FEM matrix...");fflush(stdout);}

  double tic0 = MPI_Wtime();
  
  switch(mesh.elementType){
  case TRIANGLES:
    BuildOperatorMatrixContinuousTri2D(AL); break;
  case QUADRILATERALS:
    {
      if(mesh.dim==2)
	BuildOperatorMatrixContinuousQuad2D(AL);
      else
	BuildOperatorMatrixContinuousQuad3D(AL);

      break;
    }
  case TETRAHEDRA:
    BuildOperatorMatrixContinuousTet3D(AL); break;
  case HEXAHEDRA:
    BuildOperatorMatrixContinuousHex3D(AL); break;
  }

  double tic1 = MPI_Wtime();
  printf("Local matrices took %g secs on HOST\n", tic1-tic0);
  
  compressMatrix(mesh, ogsMasked, AL, nnzLocal, A);
}

void elliptic_t::BuildOperatorMatrixContinuousTri2D(parAlmond::parCOO::nonZero_t *AL) {


  // Build non-zeros of stiffness matrix (unassembled)
  dfloat *Srr = mesh.Srr;
  dfloat *Srs = mesh.Srs;
  dfloat *Sss = mesh.Sss;
  dfloat *MM  = mesh.MM ;

  //Build unassembed non-zeros
  for (dlong e=0;e<mesh.Nelements;e++) {
    dfloat Grr = mesh.ggeo[e*mesh.Nggeo + G00ID];
    dfloat Grs = mesh.ggeo[e*mesh.Nggeo + G01ID];
    dfloat Gss = mesh.ggeo[e*mesh.Nggeo + G11ID];
    dfloat J   = mesh.ggeo[e*mesh.Nggeo + GWJID];

    for (int n=0;n<mesh.Np;n++) {
      if (maskedGlobalNumbering[e*mesh.Np + n]<0) continue; //skip masked nodes
      for (int m=0;m<mesh.Np;m++) {
        if (maskedGlobalNumbering[e*mesh.Np + m]<0) continue; //skip masked nodes

        dfloat val = 0.;

        val += Grr*Srr[m+n*mesh.Np];
        val += Grs*Srs[m+n*mesh.Np];
        val += Gss*Sss[m+n*mesh.Np];
        val += J*lambda*MM[m+n*mesh.Np];

	// pack non-zero
	dlong cnt = e*mesh.Np*mesh.Np + n*mesh.Np + m;
	AL[cnt].val = val;
	AL[cnt].row = maskedGlobalNumbering[e*mesh.Np + n];
	AL[cnt].col = maskedGlobalNumbering[e*mesh.Np + m];
      }
    }
  }

}


void elliptic_t::BuildOperatorMatrixContinuousQuad3D(parAlmond::parCOO::nonZero_t *AL) {

  //Build unassembed non-zeros
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int ny=0;ny<mesh.Nq;ny++) {
      for (int nx=0;nx<mesh.Nq;nx++) {
        if (maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq]<0) continue; //skip masked nodes
        for (int my=0;my<mesh.Nq;my++) {
          for (int mx=0;mx<mesh.Nq;mx++) {
            if (maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq]<0) continue; //skip masked nodes

            int id;
            dfloat val = 0.;

            if (ny==my) {
              for (int k=0;k<mesh.Nq;k++) {
                id = k+ny*mesh.Nq;
                dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G00ID*mesh.Np];

                val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
              }
            }

            id = mx+ny*mesh.Nq;
            dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
            val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];

            id = nx+my*mesh.Nq;
            dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
            val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];

            if (nx==mx) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+k*mesh.Nq;
                dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G11ID*mesh.Np];

                val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
              }
            }


            if ((nx==mx)&&(ny==my)) {
              id = nx + ny*mesh.Nq;

              // dfloat Gtt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G22ID*mesh.Np];
              // val += Gtt;

              dfloat JW = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + GWJID*mesh.Np];
              val += JW*lambda;
            }

	    // pack non-zero
	    dlong cnt = e*mesh.Np*mesh.Np + (nx+ny*mesh.Nq)*mesh.Np + mx+my*mesh.Nq;
	    AL[cnt].val = val;
	    AL[cnt].row = maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq];
	    AL[cnt].col = maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq];
          }
        }
      }
    }
  }

}


void elliptic_t::BuildOperatorMatrixContinuousQuad2D(parAlmond::parCOO::nonZero_t *AL) {

  //Build unassembed non-zeros
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int ny=0;ny<mesh.Nq;ny++) {
      for (int nx=0;nx<mesh.Nq;nx++) {
        if (maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq]<0) continue; //skip masked nodes
        for (int my=0;my<mesh.Nq;my++) {
          for (int mx=0;mx<mesh.Nq;mx++) {
            if (maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq]<0) continue; //skip masked nodes

            int id;
            dfloat val = 0.;

            if (ny==my) {
              for (int k=0;k<mesh.Nq;k++) {
                id = k+ny*mesh.Nq;
                dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G00ID*mesh.Np];

                val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
              }
            }

            id = mx+ny*mesh.Nq;
            dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
            val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];


            id = nx+my*mesh.Nq;
            dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
            val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];

            if (nx==mx) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+k*mesh.Nq;
                dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G11ID*mesh.Np];

                val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
              }
            }

            if ((nx==mx)&&(ny==my)) {
              id = nx + ny*mesh.Nq;
              dfloat JW = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + GWJID*mesh.Np];
              val += JW*lambda;
            }

	    // pack non-zero
	    dlong cnt = e*mesh.Np*mesh.Np + (nx+ny*mesh.Nq)*mesh.Np + mx+my*mesh.Nq;
	    AL[cnt].val = val;
	    AL[cnt].row = maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq];
	    AL[cnt].col = maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq];
          }
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorMatrixContinuousTet3D(parAlmond::parCOO::nonZero_t *AL) {

  for (dlong e=0;e<mesh.Nelements;e++) {

    dfloat Grr = mesh.ggeo[e*mesh.Nggeo + G00ID];
    dfloat Grs = mesh.ggeo[e*mesh.Nggeo + G01ID];
    dfloat Grt = mesh.ggeo[e*mesh.Nggeo + G02ID];
    dfloat Gss = mesh.ggeo[e*mesh.Nggeo + G11ID];
    dfloat Gst = mesh.ggeo[e*mesh.Nggeo + G12ID];
    dfloat Gtt = mesh.ggeo[e*mesh.Nggeo + G22ID];
    dfloat J   = mesh.ggeo[e*mesh.Nggeo + GWJID];

    for (int n=0;n<mesh.Np;n++) {
      if (maskedGlobalNumbering[e*mesh.Np + n]<0) continue; //skip masked nodes
      for (int m=0;m<mesh.Np;m++) {
        if (maskedGlobalNumbering[e*mesh.Np + m]<0) continue; //skip masked nodes
        dfloat val = 0.;

        val += Grr*mesh.Srr[m+n*mesh.Np];
        val += Grs*mesh.Srs[m+n*mesh.Np];
        val += Grt*mesh.Srt[m+n*mesh.Np];
        val += Gss*mesh.Sss[m+n*mesh.Np];
        val += Gst*mesh.Sst[m+n*mesh.Np];
        val += Gtt*mesh.Stt[m+n*mesh.Np];
        val += J*lambda*mesh.MM[m+n*mesh.Np];

	// pack non-zero
	dlong cnt = e*mesh.Np*mesh.Np + n*mesh.Np + m;
	AL[cnt].val = val;
	AL[cnt].row = maskedGlobalNumbering[e*mesh.Np + n];
	AL[cnt].col = maskedGlobalNumbering[e*mesh.Np + m];
      }
    }
  }
}

void elliptic_t::BuildOperatorMatrixContinuousHex3D(parAlmond::parCOO::nonZero_t *AL) {

  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int nz=0;nz<mesh.Nq;nz++) {
      for (int ny=0;ny<mesh.Nq;ny++) {
	for (int nx=0;nx<mesh.Nq;nx++) {
	  int idn = nx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
	  if (maskedGlobalNumbering[e*mesh.Np + idn]<0) continue; //skip masked nodes
	
	  for (int mz=0;mz<mesh.Nq;mz++) {
	    for (int my=0;my<mesh.Nq;my++) {
	      for (int mx=0;mx<mesh.Nq;mx++) {
		int idm = mx+my*mesh.Nq+mz*mesh.Nq*mesh.Nq;
		if (maskedGlobalNumbering[e*mesh.Np + idm]<0) continue; //skip masked nodes
	      
		int id;
		dfloat val = 0.;
	      
		if ((ny==my)&&(nz==mz)) {
		  for (int k=0;k<mesh.Nq;k++) {
		    id = k+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		    dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G00ID*mesh.Np];
		  
		    val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
		  }
		}
	      
		if (nz==mz) {
		  id = mx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		  dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
		  val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];
		
		  id = nx+my*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		  dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G01ID*mesh.Np];
		  val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];
		}
	      
		if (ny==my) {
		  id = mx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		  dfloat Grt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G02ID*mesh.Np];
		  val += Grt*mesh.D[nx+mx*mesh.Nq]*mesh.D[mz+nz*mesh.Nq];
		
		  id = nx+ny*mesh.Nq+mz*mesh.Nq*mesh.Nq;
		  dfloat Gst = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G02ID*mesh.Np];
		  val += Gst*mesh.D[mx+nx*mesh.Nq]*mesh.D[nz+mz*mesh.Nq];
		}
	      
		if ((nx==mx)&&(nz==mz)) {
		  for (int k=0;k<mesh.Nq;k++) {
		    id = nx+k*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		    dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G11ID*mesh.Np];
		  
		    val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
		  }
		}
	      
		if (nx==mx) {
		  id = nx+my*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		  dfloat Gst = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G12ID*mesh.Np];
		  val += Gst*mesh.D[ny+my*mesh.Nq]*mesh.D[mz+nz*mesh.Nq];
		
		  id = nx+ny*mesh.Nq+mz*mesh.Nq*mesh.Nq;
		  dfloat Gts = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G12ID*mesh.Np];
		  val += Gts*mesh.D[my+ny*mesh.Nq]*mesh.D[nz+mz*mesh.Nq];
		}

		if ((nx==mx)&&(ny==my)) {
		  for (int k=0;k<mesh.Nq;k++) {
		    id = nx+ny*mesh.Nq+k*mesh.Nq*mesh.Nq;
		    dfloat Gtt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + G22ID*mesh.Np];

		    val += Gtt*mesh.D[nz+k*mesh.Nq]*mesh.D[mz+k*mesh.Nq];
		  }
		}

		if ((nx==mx)&&(ny==my)&&(nz==mz)) {
		  id = nx + ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
		  dfloat JW = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + GWJID*mesh.Np];
		  val += JW*lambda;
		}

		dlong cnt = e*mesh.Np*mesh.Np +
		  idn*mesh.Np +
		  idm;

		AL[cnt].val = val;
		AL[cnt].row = maskedGlobalNumbering[e*mesh.Np + idn];
		AL[cnt].col = maskedGlobalNumbering[e*mesh.Np + idm];
	      }
	    }
	  }
	}
      }
    }
  }

}
