#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "mesh3D.h"

// interpolate data to plot nodes and save to file (one per process
void meshPlotVTU3DP(mesh3D *mesh, char *fileNameBase, iint fld){
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *fp;
  char fileName[BUFSIZ];
  //sprintf(fileName, "%s_%04d.vtu", fileNameBase, rank);
  strcpy(fileName,fileNameBase);
  
  iint NMax = mesh->NMax;

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	  mesh->Nelements*mesh->plotNp[NMax], 
	  mesh->Nelements*mesh->plotNelements[NMax]);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(iint e=0;e<mesh->Nelements;++e){
    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    
    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    for(iint n=0;n<mesh->plotNp[NMax];++n){
      /* (r,s,t) coordinates of plot nodes*/
      dfloat rn = mesh->plotR[NMax][n]; 
      dfloat sn = mesh->plotS[NMax][n];
      dfloat tn = mesh->plotT[NMax][n];

      /* physical coordinate of interpolation node */
      dfloat plotxn = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
      dfloat plotyn = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
      dfloat plotzn = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
      
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" Format=\"ascii\">\n");
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp[NMax];++n){
      dfloat plotpn = 0;
      for(iint m=0;m<mesh->NpMax;++m){
        dfloat pm = mesh->q[fld + mesh->Nfields*(m+e*mesh->NpMax)];
        plotpn += mesh->plotInterp[NMax][n*mesh->NpMax+m]*pm;
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }

  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements[NMax];++n){
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
	fprintf(fp, "%d ", e*mesh->plotNp[NMax] + mesh->plotEToV[NMax][n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements[NMax];++n){
      cnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements[NMax];++n){
      fprintf(fp, "10\n"); // TET code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
