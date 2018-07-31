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

// holmes/solvers/bns/src/bnsWriteIsoGmsh.c
// interpolate data to plot nodes and save to file (one per process)
// 2018/06/02
//---------------------------------------------------------
#include "bns.h"

void bnsIsoPlotGmsh(
  bns_t *bns, 
  int Ntris1, 
  char *fname,
  int tstep,        // simulation time-step
  int N_offset,     // gmsh mpi node offset
  int E_offset,     // gmsh mpi element offset
  int plotnum,      // plot counter
  double plottime,  // simulation time
  bool bBinary)
{

  // Export isosurf data interpolated onto plot nodes.
  // Each proc exports its part of the solution in gmsh format, 
  // adjusting node/elem nums to allow gmsh to merge all parts.

  mesh_t *mesh = bns->mesh;
  // const int bns->procid = 0;
  // const int bns->procid = mesh->m_bns->procid;

  FILE *fp = fopen(fname,bBinary ? "wb" : "w");
  if (!fp) {
    printf("[proc:%02d] unable to open output file %s\n",bns->procid,fname);
    return;
  }

  // for isosurfaces, geometry/topology changes for each plot,
  // so we need to export mesh coord data for every plot.

  std::vector<double>& refN = bns->iso_nodes;
  std::vector<int>& refT = bns->iso_tris;

  int Nnodes = refN.size()/(bns->dim + bns->isoNfields);
  int Ntris2 = refT.size()/(bns->dim);


#if (1)
  printf("\n[proc:%02d] Report from bnsWriteIsoGmsh()\n"
             "  - total tris  to export: %6d\n" 
             "  - elem offset for  proc: %6d\n" 
             "  - node offset for  proc: %6d\n", bns->procid, 
             Ntris1, E_offset, N_offset);
#endif

  //-------------------------------------
  // export in gmsh 2.2 format
  //-------------------------------------

  // size of data on this platform
  size_t szD = sizeof(double), szI = sizeof(int);

  //-------------------------------------
  // 1. Write the gmsh header details
  //-------------------------------------
  fprintf(fp, "$MeshFormat\n");             // [Section] start
  if (bBinary) {
    fprintf(fp, "2.2  1  %d\n", (int)szD);  // binary
    int one=1; fwrite(&one, szI, 1, fp);    // test binary "endian"
    fprintf(fp, "\n");                      // add eol for binary
  } else {
    fprintf(fp, "2.2  0  %d\n", (int)szD);  // ascii
  }
  fprintf(fp, "$EndMeshFormat \n");         // [Section] stop


  // Start of geometry/topology data

  //-------------------------------------
  // 2. Write the vertex data
  //-------------------------------------
  fprintf(fp,"$Nodes\n");             // [Section] start
  fprintf(fp,"%d\n", Nnodes);         // total nodes

  //                 3             1
  int noff = (mesh->dim+bns->isoNfields);
  int id1 = N_offset;  // node offset for this partition

  if (bBinary) {
    double data[3];
    for(dlong n=0; n<Nnodes; ++n){

      int idn = n*noff;         // load coords for node n
      data[0] = refN[idn+0];
      data[1] = refN[idn+1];
      data[2] = refN[idn+2];

      ++id1;                    // make node id 1-based for gmsh
      fwrite(&id1, szI,1,fp);
      fwrite(data, szD,3,fp);
    }
    fprintf(fp,"\n");           // add eol for binary
  }
  else {
    double pxn,pyn,pzn;
    for (dlong n=0; n<Nnodes; ++n) {

      int idn = n*noff;         // load coords for node n
      pxn = refN[idn+0]; 
      pyn = refN[idn+1];
      pzn = refN[idn+2];

      // make node id 1-based for gmsh
      fprintf(fp,"%7d  %14.6e %14.6e %14.6e\n", ++id1, pxn,pyn,pzn);
    }
  }
  fprintf(fp,"$EndNodes\n");    // [Section] stop


  //-------------------------------------
  // 3. Write element connectivity
  //-------------------------------------
  fprintf(fp,"$Elements\n");    // [Section] start 
  fprintf(fp,"%d\n",Ntris1);    // total elements

  id1 = E_offset;               // element offset for this partition
  int nodesk = N_offset;        // mpi proc offset for gmsh

  //------------------------------------------------------
  //                         0   1 2 3 4         5 6 7 [8]
  int Ntags=4; int eData[9]={0, 99,2,1,bns->procid+1, 0,0,0, 0 };
  //------------------------------------------------------

  if (bBinary) {
    // write header for group of tri elements
    int header[3] = { 2,Ntris1,Ntags }; fwrite(&header,szI,3,fp);
  }

  // for each tri element ...
  int n1,n2,n3;
  for (dlong e=0; e<Ntris1; ++e){

    n1 = nodesk + refT[e*3 +0] +1;    // make ids 1-based for gmsh
    n2 = nodesk + refT[e*3 +1] +1;
    n3 = nodesk + refT[e*3 +2] +1;

    if (bBinary) {
      eData[0] = ++id1;               // increment/load element id
      eData[5] = n1;                  // load 3 ids for next tri
      eData[6] = n2;
      eData[7] = n3;                  //       id   4 tags   3 nodes
      fwrite(eData,szI,8,fp);         // e.g. [id, 99,2,1,3  a,b,c]
    }
    else {        // [tp:2] ==> tri
      // write   [id][tp][Nt,t1,2,3,4] + [3] ids for this tri
      fprintf(fp,"%7d 2    4 99 2 1 %d  %6d %6d %6d\n",
                  ++id1,      bns->procid+1, n1, n2, n3);
    }
  }

  if (bBinary) { fprintf(fp,"\n"); }    // add eol for binary
  fprintf(fp,"$EndElements\n");         // [Section] stop

  // End of geometry/topology data


  //-------------------------------------
  // 4. Write the (scalar) node data
  //-------------------------------------

  // For each field, write NODE data for each node
  // int num_fields = 1;

  {
    int numComp = 1;                            // 1-component: scalar field
    for (int fld=1; fld<=bns->isoNfields; ++fld) {   // [Section] start 
      fprintf(fp, "$NodeData\n");               // scalar data for next field
      fprintf(fp, "1\n");                       // one string tag
      fprintf(fp, "\"%s\"\n", "iso");           // displayName (MUST double quote)
      fprintf(fp, "1\n%0.16g\n", plottime);     // one real tag: time 
      fprintf(fp, "4\n%d\n%d\n%d\n%d\n",        // 4 int tags:
                   plotnum, numComp, Nnodes,    // step,Ncomp,Nnode
                   bns->procid+1);                   // partition (1-based)

      int id = N_offset;   // node offset for this partition
      if (bBinary) {
        double val=0.0;

        for (dlong n=0; n<Nnodes; ++n) {
          ++id; val = refN[n*noff + fld];   // incr id, get value
          fwrite(&id,szI,1,fp);           // write id
          fwrite(&val,szD,1,fp);          // write value
        }
        fprintf(fp, "\n");                // add eol for binary
      }
      else {
        for (dlong n=0; n<Nnodes; ++n) {
          fprintf(fp, "%d  %14.6e \n", ++id, refN[n*noff + fld]);
        }
      }
      fprintf(fp, "$EndNodeData\n");      // [Section] stop
    }
  }

  fclose(fp);
}

