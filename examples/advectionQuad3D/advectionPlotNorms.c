#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "advectionQuad3D.h"

// interpolate data to plot nodes and save to file (one per process
void advectionPlotNorms(mesh_t *mesh, char *fileNameBase, int tstep,dfloat *q){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  FILE *fp;
  char fileName[BUFSIZ];
  sprintf(fileName, "%s_%04d_%04d.vtu", fileNameBase, rank, tstep);
  
  printf("FILE = %s\n", fileName);

  fp = fopen(fileName, "w");
  
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	  mesh->Nelements*mesh->Np,
	  mesh->Nelements*mesh->Np);

  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat plotxn = mesh->x[n+e*mesh->Np];
      dfloat plotyn = mesh->y[n+e*mesh->Np];
      dfloat plotzn = mesh->z[n+e*mesh->Np];
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"mrab_levels\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat plotpn = mesh->MRABlevel[e];
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }

  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"q_1\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      dfloat plotpn = q[n+e*mesh->Np];
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"jacobian\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n) {
      int id = mesh->Nvgeo*mesh->Np*e + JWID*mesh->Np + n;
      dfloat plotpn = mesh->vgeo[id];
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "    </PointData>\n");
      
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<mesh->Nelements;++e){
	fprintf(fp, "       ");
	fprintf(fp, "%d %d %d %d\n", e*mesh->Np + 0,e*mesh->Np + mesh->Nq - 1,(e+1)*mesh->Np-1,(e+1)*mesh->Np - mesh->Nq);
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
      cnt += 4;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
      fprintf(fp, "9\n");
  }

  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}
