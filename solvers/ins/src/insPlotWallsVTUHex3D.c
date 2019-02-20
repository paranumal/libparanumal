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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "ins.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotWallsVTUHex3D(ins_t *ins, char *fileNameBase){

  mesh_t *mesh = ins->mesh;

  int rank;
  rank = mesh->rank;

  // count walls
  hlong Nwalls = 0;
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToB[e*mesh->Nfaces+f]==1) { // need to introduce defines
	++Nwalls;
      }
    }
  }

  FILE *fp;
  char fileName[BUFSIZ];
  sprintf(fileName, "%s_%04d.vtu", fileNameBase, rank);

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");

  printf("N = %d, Nwalls = " hlongFormat ", Nel = %d\n",
	 mesh->Nq-1, Nwalls, mesh->Nelements);

  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"" hlongFormat "\">\n",
          mesh->Nelements*mesh->Np,
          Nwalls*(mesh->Nq-1)*(mesh->Nq-1));

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  // compute plot node coordinates on the fly
  for(hlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      hlong id = n + e*mesh->Np;
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n",
	      mesh->x[id],
	      mesh->y[id],
	      mesh->z[id]);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToB[e*mesh->Nfaces+f]==1) { // need to introduce defines
	hlong b = e*mesh->Np;
	for(int j=0;j<mesh->Nq-1;++j){
	  for(int i=0;i<mesh->Nq-1;++i){
	    int v1 = mesh->faceNodes[f*mesh->Nfp + j*mesh->Nq + i];
	    int v2 = mesh->faceNodes[f*mesh->Nfp + j*mesh->Nq + i+1];
	    int v3 = mesh->faceNodes[f*mesh->Nfp + (j+1)*mesh->Nq + i+1];
	    int v4 = mesh->faceNodes[f*mesh->Nfp + (j+1)*mesh->Nq + i];

	    fprintf(fp,
		    hlongFormat" "
		    hlongFormat" "
		    hlongFormat" "
		    hlongFormat"\n ",
		    b + v1, b + v2, b + v3, b + v4);
	  }
	}
      }
    }
  }

  fprintf(fp, "        </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong w=0;w<Nwalls*(mesh->Nq-1)*(mesh->Nq-1);++w){
    cnt += 4;
    fprintf(fp, "       ");
    fprintf(fp, dlongFormat"\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong w=0;w<Nwalls*(mesh->Nq-1)*(mesh->Nq-1);++w){
    fprintf(fp, "9\n"); // quad code ?
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
