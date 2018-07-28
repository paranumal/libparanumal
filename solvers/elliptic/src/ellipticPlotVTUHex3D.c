#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "mesh3D.h"

// interpolate data to plot nodes and save to file (one per process
extern "C"
{
void ellipticPlotVTUHex3D(mesh3D *mesh, char *fileNameBase, int fld);
}

void ellipticPlotVTUHex3D(mesh3D *mesh, char *fileNameBase, int fld){
  
  int rank;
  rank = mesh->rank;

  FILE *fp;
  char fileName[BUFSIZ];
  sprintf(fileName, "%s_%04d.vtu", fileNameBase, rank);
  // strcpy(fileName,fileNameBase);

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");

  int Eloc = (mesh->Nq-1)*(mesh->Nq-1)*(mesh->Nq-1);
  printf("N = %d, Eloc = %d, Nel = %d\n",
	 mesh->Nq-1, Eloc, mesh->Nelements);

  fprintf(fp, "    <Piece NumberOfPoints=\""dlongFormat"\" NumberOfCells=\""dlongFormat"\">\n", 
          mesh->Nelements*mesh->Np, 
          mesh->Nelements*Eloc);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dlong id = n + e*mesh->Np;
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n",
	      mesh->x[id],
	      mesh->y[id],
	      mesh->z[id]);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"pressure\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat qn = mesh->q[n+fld*mesh->Np + e*mesh->Nfields*mesh->Np];
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", qn);
    }
  }

  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int k=0;k<mesh->Nq-1;++k){
      for(int j=0;j<mesh->Nq-1;++j){
	for(int i=0;i<mesh->Nq-1;++i){
	  int b = e*mesh->Np + i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;
	  fprintf(fp, 
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat" "
		  dlongFormat"\n ",
		  b,
		  b+1,
		  b+mesh->Nq+1,
		  b+mesh->Nq,
		  b + mesh->Nq*mesh->Nq,
		  b+1 + mesh->Nq*mesh->Nq,
		  b+mesh->Nq+1 + mesh->Nq*mesh->Nq,
		  b+mesh->Nq+ mesh->Nq*mesh->Nq);
	}
      }
    }
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<Eloc;++n){
      cnt += 8;
      fprintf(fp, "       ");
      fprintf(fp, dlongFormat"\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<Eloc;++n){
      fprintf(fp, "12\n"); // HEX code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}
