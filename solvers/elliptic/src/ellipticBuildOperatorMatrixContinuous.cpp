/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

void elliptic_t::BuildOperatorMatrixContinuous(parAlmond::parCOO& A) {

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
    BuildOperatorMatrixContinuousTri2D(A); break;
  case Mesh::QUADRILATERALS:
  {
    if(mesh.dim==2)
      BuildOperatorMatrixContinuousQuad2D(A);
    else
      BuildOperatorMatrixContinuousQuad3D(A);

    break;
  }
  case Mesh::TETRAHEDRA:
    BuildOperatorMatrixContinuousTet3D(A); break;
  case Mesh::HEXAHEDRA:
    BuildOperatorMatrixContinuousHex3D(A); break;
  }
}

void elliptic_t::BuildOperatorMatrixContinuousTri2D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1, 0);
  A.globalColStarts.malloc(mesh.size+1, 0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;

  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  memory<dfloat> Srr = mesh.Srr;
  memory<dfloat> Srs = mesh.Srs;
  memory<dfloat> Sss = mesh.Sss;
  memory<dfloat> MM  = mesh.MM ;

  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

  //Build unassembed non-zeros
  dlong cnt =0;
  for (dlong e=0;e<mesh.Nelements;e++) {
    dfloat Grr = mesh.ggeo[e*mesh.Nggeo + mesh.G00ID];
    dfloat Grs = mesh.ggeo[e*mesh.Nggeo + mesh.G01ID];
    dfloat Gss = mesh.ggeo[e*mesh.Nggeo + mesh.G11ID];
    dfloat J   = mesh.wJ[e];

    for (int n=0;n<mesh.Np;n++) {
      if (maskedGlobalNumbering[e*mesh.Np + n]<0) continue; //skip masked nodes
      for (int m=0;m<mesh.Np;m++) {
        if (maskedGlobalNumbering[e*mesh.Np + m]<0) continue; //skip masked nodes

        dfloat val = 0.;

        val += Grr*Srr[m+n*mesh.Np];
        val += Grs*Srs[m+n*mesh.Np];
        val += Gss*Sss[m+n*mesh.Np];
        val += J*lambda*MM[m+n*mesh.Np];

        dfloat nonZeroThreshold = 1e-7;
        if (fabs(val)>nonZeroThreshold) {
          // pack non-zero
          sendNonZeros[cnt].val = val;
          sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + n];
          sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + m];
          cnt++;
        }
      }
    }
  }

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
       [](const parAlmond::parCOO::nonZero_t& a,
          const parAlmond::parCOO::nonZero_t& b) {
         if (a.row < b.row) return true;
         if (a.row > b.row) return false;

         return a.col < b.col;
        });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
       [](const parAlmond::parCOO::nonZero_t& a,
          const parAlmond::parCOO::nonZero_t& b) {
         if (a.row < b.row) return true;
         if (a.row > b.row) return false;

         return a.col < b.col;
       });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

  if(Comm::World().rank()==0) printf("done.\n");
}


void elliptic_t::BuildOperatorMatrixContinuousQuad3D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1, 0);
  A.globalColStarts.malloc(mesh.size+1, 0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;
  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

#if 0
  hlong NTf = mesh.Nelements*mesh.Np * mesh.Nelements*mesh.Np ;
  dfloat *Af = (dfloat *)calloc(NTf, sizeof(dfloat));
#endif

  //Build unassembed non-zeros
  dlong cnt =0;
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
                dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G00ID*mesh.Np];

                val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
              }
            }

            id = mx+ny*mesh.Nq;
            dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
            val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];

            id = nx+my*mesh.Nq;
            dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
            val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];


            // id = mx+ny*mesh.Nq;
            // dfloat Grt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G02ID*mesh.Np];
            // val += Grt*mesh.D[nx+mx*mesh.Nq];

            // id = nx+my*mesh.Nq;
            // dfloat Gtr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G02ID*mesh.Np];
            // val += Gtr*mesh.D[mx+nx*mesh.Nq];


            if (nx==mx) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+k*mesh.Nq;
                dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G11ID*mesh.Np];

                val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
              }
            }

            // double check following two: AK
            // id = nx+my*mesh.Nq;
            // dfloat Gst = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G12ID*mesh.Np];
            // val += Gst*mesh.D[ny+my*mesh.Nq];

            // id = mx+ny*mesh.Nq;
            // dfloat Gts = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G12ID*mesh.Np];
            // val += Gts*mesh.D[my+ny*mesh.Nq];


            if ((nx==mx)&&(ny==my)) {
              id = nx + ny*mesh.Nq;

              // dfloat Gtt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G22ID*mesh.Np];
              // val += Gtt;

              dfloat JW = mesh.wJ[e*mesh.Np + id];
              val += JW*lambda;
            }

#if 0
            const hlong rowid = e*mesh.Np + nx + ny*mesh.Nq;
            const hlong colid = e*mesh.Np + mx + my*mesh.Nq;

            Af[rowid*mesh.Nelements*mesh.Np + colid] = val;
#endif

            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val)>nonZeroThreshold) {
              // pack non-zero
              sendNonZeros[cnt].val = val;
              sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq];
              sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq];
              cnt++;
            }
          }
        }
      }
    }
  }

#if 0
 // Write matlab dat for postprocess
  char fname[BUFSIZ];
  sprintf(fname, "Ax.dat");
  FILE *fp;
  fp = fopen(fname, "w");

  for(hlong row = 0; row<(mesh.Nelements*mesh.Np); row++){
    for(hlong col = 0; col<(mesh.Nelements*mesh.Np); col++){
      dfloat val = Af[row*mesh.Nelements*mesh.Np + col];
      fprintf(fp,"%.8e ", val);
    }
    fprintf(fp,"\n");
  }

 fclose(fp);
#endif

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
        [](const parAlmond::parCOO::nonZero_t& a,
           const parAlmond::parCOO::nonZero_t& b) {
          if (a.row < b.row) return true;
          if (a.row > b.row) return false;

          return a.col < b.col;
        });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
        [](const parAlmond::parCOO::nonZero_t& a,
           const parAlmond::parCOO::nonZero_t& b) {
          if (a.row < b.row) return true;
          if (a.row > b.row) return false;

          return a.col < b.col;
        });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

#if 0
  // Write matlab dat for postprocess
  char fname[BUFSIZ];
  sprintf(fname, "Ax.dat");
  FILE *fp;
  fp = fopen(fname, "w");

  for(dlong n=1;n<A.nnz;++n){
      fprintf(fp,"%d %d %.8e\n", (*A)[n].row+1, (*A)[n].col+1, (*A)[n].val);
  }

 fclose(fp);
#endif

  if(Comm::World().rank()==0) printf("done.\n");
}


void elliptic_t::BuildOperatorMatrixContinuousQuad2D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;
  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

  //Build unassembed non-zeros
  dlong cnt =0;
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
                dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G00ID*mesh.Np];

                val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
              }
            }

            id = mx+ny*mesh.Nq;
            dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
            val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];


            id = nx+my*mesh.Nq;
            dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
            val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];

            if (nx==mx) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+k*mesh.Nq;
                dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G11ID*mesh.Np];

                val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
              }
            }

            if ((nx==mx)&&(ny==my)) {
              id = nx + ny*mesh.Nq;
              dfloat JW = mesh.wJ[e*mesh.Np + id];
              val += JW*lambda;
            }

            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val)>nonZeroThreshold) {
              // pack non-zero
              sendNonZeros[cnt].val = val;
              sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + nx+ny*mesh.Nq];
              sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + mx+my*mesh.Nq];
              cnt++;
            }
          }
        }
      }
    }
  }

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

#if 0
  // Write matlab dat for postprocess
  char fname[BUFSIZ];
  sprintf(fname, "Ax.dat");
  FILE *fp;
  fp = fopen(fname, "w");

  for(dlong n=1;n<A.nnz;++n){
      fprintf(fp, hlongFormat " " hlongFormat " %.8e\n", (*A)[n].row+1, (*A)[n].col+1, (*A)[n].val);
  }

 fclose(fp);
#endif

  if(Comm::World().rank()==0) printf("done.\n");
}

void elliptic_t::BuildOperatorMatrixContinuousTet3D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;

  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  //Build unassembed non-zeros
  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

  dlong cnt =0;
  //#pragma omp parallel for
  for (dlong e=0;e<mesh.Nelements;e++) {

    dfloat Grr = mesh.ggeo[e*mesh.Nggeo + mesh.G00ID];
    dfloat Grs = mesh.ggeo[e*mesh.Nggeo + mesh.G01ID];
    dfloat Grt = mesh.ggeo[e*mesh.Nggeo + mesh.G02ID];
    dfloat Gss = mesh.ggeo[e*mesh.Nggeo + mesh.G11ID];
    dfloat Gst = mesh.ggeo[e*mesh.Nggeo + mesh.G12ID];
    dfloat Gtt = mesh.ggeo[e*mesh.Nggeo + mesh.G22ID];
    dfloat J   = mesh.wJ[e];

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

        dfloat nonZeroThreshold = 1e-7;
        if (fabs(val)>nonZeroThreshold) {
          //#pragma omp critical
          {
            // pack non-zero
            sendNonZeros[cnt].val = val;
            sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + n];
            sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + m];
            cnt++;
          }
        }
      }
    }
  }

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

  if(Comm::World().rank()==0) printf("done.\n");
}

void elliptic_t::BuildOperatorMatrixContinuousHex3D(parAlmond::parCOO& A) {

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = ogsMasked.Ngather;

  // every gathered degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Ngather, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  // 2. Build non-zeros of stiffness matrix (unassembled)
  dlong nnzLocal = mesh.Np*mesh.Np*mesh.Nelements;
  memory<parAlmond::parCOO::nonZero_t> sendNonZeros(nnzLocal);
  memory<int> AsendCounts (mesh.size, 0);
  memory<int> ArecvCounts (mesh.size);
  memory<int> AsendOffsets(mesh.size+1);
  memory<int> ArecvOffsets(mesh.size+1);

  if(Comm::World().rank()==0) {printf("Building full FEM matrix...");fflush(stdout);}

  dlong cnt =0;
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
                dfloat Grr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G00ID*mesh.Np];

                val += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[mx+k*mesh.Nq];
              }
            }

            if (nz==mz) {
              id = mx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
              dfloat Grs = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
              val += Grs*mesh.D[nx+mx*mesh.Nq]*mesh.D[my+ny*mesh.Nq];

              id = nx+my*mesh.Nq+nz*mesh.Nq*mesh.Nq;
              dfloat Gsr = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
              val += Gsr*mesh.D[mx+nx*mesh.Nq]*mesh.D[ny+my*mesh.Nq];
            }

            if (ny==my) {
              id = mx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
              dfloat Grt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G02ID*mesh.Np];
              val += Grt*mesh.D[nx+mx*mesh.Nq]*mesh.D[mz+nz*mesh.Nq];

              id = nx+ny*mesh.Nq+mz*mesh.Nq*mesh.Nq;
              dfloat Gst = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G02ID*mesh.Np];
              val += Gst*mesh.D[mx+nx*mesh.Nq]*mesh.D[nz+mz*mesh.Nq];
            }

            if ((nx==mx)&&(nz==mz)) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+k*mesh.Nq+nz*mesh.Nq*mesh.Nq;
                dfloat Gss = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G11ID*mesh.Np];

                val += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[my+k*mesh.Nq];
              }
            }

            if (nx==mx) {
              id = nx+my*mesh.Nq+nz*mesh.Nq*mesh.Nq;
              dfloat Gst = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G12ID*mesh.Np];
              val += Gst*mesh.D[ny+my*mesh.Nq]*mesh.D[mz+nz*mesh.Nq];

              id = nx+ny*mesh.Nq+mz*mesh.Nq*mesh.Nq;
              dfloat Gts = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G12ID*mesh.Np];
              val += Gts*mesh.D[my+ny*mesh.Nq]*mesh.D[nz+mz*mesh.Nq];
            }

            if ((nx==mx)&&(ny==my)) {
              for (int k=0;k<mesh.Nq;k++) {
                id = nx+ny*mesh.Nq+k*mesh.Nq*mesh.Nq;
                dfloat Gtt = mesh.ggeo[e*mesh.Np*mesh.Nggeo + id + mesh.G22ID*mesh.Np];

                val += Gtt*mesh.D[nz+k*mesh.Nq]*mesh.D[mz+k*mesh.Nq];
              }
            }

            if ((nx==mx)&&(ny==my)&&(nz==mz)) {
              id = nx + ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
              dfloat JW = mesh.wJ[e*mesh.Np + id];
              val += JW*lambda;
            }

            // pack non-zero
            dfloat nonZeroThreshold = 1e-7;
            if (fabs(val) >= nonZeroThreshold) {
              sendNonZeros[cnt].val = val;
              sendNonZeros[cnt].row = maskedGlobalNumbering[e*mesh.Np + idn];
              sendNonZeros[cnt].col = maskedGlobalNumbering[e*mesh.Np + idm];
              cnt++;
            }
        }
        }
        }
      }
      }
      }
  }

  // sort by row ordering
  sort(sendNonZeros.ptr(), sendNonZeros.ptr()+cnt,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // count how many non-zeros to send to each process
  int rr=0;
  for(dlong n=0;n<cnt;++n) {
    const hlong id = sendNonZeros[n].row;
    while(id>=A.globalRowStarts[rr+1]) rr++;
    AsendCounts[rr]++;
  }

  // find how many nodes to expect (should use sparse version)
  mesh.comm.Alltoall(AsendCounts, ArecvCounts);

  // find send and recv offsets for gather
  A.nnz = 0;
  AsendOffsets[0] = 0;
  ArecvOffsets[0] = 0;
  for(int r=0;r<mesh.size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    A.nnz += ArecvCounts[r];
  }

  A.entries.malloc(A.nnz);

  // determine number to receive
  mesh.comm.Alltoallv(sendNonZeros, AsendCounts, AsendOffsets,
                      A.entries,    ArecvCounts, ArecvOffsets);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  sort(A.entries.ptr(), A.entries.ptr()+A.nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  // compress duplicates
  cnt = 0;
  for(dlong n=1;n<A.nnz;++n){
    if(A.entries[n].row == A.entries[cnt].row &&
       A.entries[n].col == A.entries[cnt].col){
       A.entries[cnt].val += A.entries[n].val;
    }
    else{
      ++cnt;
      A.entries[cnt] = A.entries[n];
    }
  }
  if (A.nnz) cnt++;
  A.nnz = cnt;

  if(Comm::World().rank()==0) printf("done.\n");
}
