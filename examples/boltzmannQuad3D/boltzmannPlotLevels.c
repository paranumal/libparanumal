#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "boltzmannQuad3D.h"

// interpolate data to plot nodes and save to file (one per process
void boltzmannPlotLevels(mesh_t *mesh, char *fileNameBase, int tstep,dfloat *q){

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
	  mesh->Nelements*mesh->plotNp, 
	  mesh->Nelements*mesh->plotNelements);

  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;
      for(iint m=0;m<mesh->Np;++m){
	plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
	plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
	plotzn += mesh->plotInterp[n*mesh->Np+m]*mesh->z[m+e*mesh->Np];
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"mrab_levels\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = mesh->MRABlevel[e];
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }

  fprintf(fp, "       </DataArray>\n");

  /*fprintf(fp, "        <DataArray type=\"Float32\" Name=\"q_value\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){

      dfloat plotpn = 0;
      for(iint m=0;m<mesh->Np;++m){
	plotpn += mesh->plotInterp[n*mesh->Np+m]*q[m+1*mesh->Np+e*mesh->Nrhs*mesh->Np*mesh->Nfields];
      }
      
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");*/

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"cube_faces\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = mesh->cubeFaceNumber[e];
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"cube_distances\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = mesh->cubeDistance[e];
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"elements\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = e%10;
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements;++n){
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
  fprintf(fp, "%d ", e*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements;++n){
      cnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements;++n){
      fprintf(fp, "5\n");
    }
  }

  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}
