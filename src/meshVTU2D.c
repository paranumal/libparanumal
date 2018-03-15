#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "mesh2D.h"

void meshVTU2D(mesh2D *mesh, char *fileName){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE *fp;
  iint *allNelements = (iint*) calloc(size, sizeof(iint));
  iint totalNelements = 0, maxNelements = 0;
  dfloat *tmpEX, *tmpEY;

  if(rank==0){
    fp = fopen(fileName, "w");
    printf("fp=%p\n for fileName=%s\n",fp, fileName);
  }

  // gather element counts to root
  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, 
		allNelements, 1, MPI_IINT, 
		MPI_COMM_WORLD);
  
  if(rank==0){

    for(int r=0;r<size;++r){
      totalNelements+=allNelements[r];
      maxNelements = mymax(maxNelements, allNelements[r]);
    }

    tmpEX = (dfloat*) calloc(maxNelements*mesh->Nverts, sizeof(dfloat));
    tmpEY = (dfloat*) calloc(maxNelements*mesh->Nverts, sizeof(dfloat));


    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", totalNelements*mesh->Nverts, totalNelements);

    // write out nodes
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

    if(rank==0){  // root writes out its coordinates
      printf("printing local verts \n");
      for(iint e=0;e<mesh->Nelements;++e){
	fprintf(fp, "        ");
	for(iint n=0;n<mesh->Nverts;++n)
	  fprintf(fp, "%g %g 0.\n", mesh->EX[e*mesh->Nverts+n], mesh->EY[e*mesh->Nverts+n]);
      }
    }
  }
  
  for(int r=1;r<size;++r){
    
    if(rank==r){
      MPI_Send(mesh->EX, mesh->Nelements*mesh->Nverts, 
	       MPI_DFLOAT, 0, 666, MPI_COMM_WORLD);
      MPI_Send(mesh->EY, mesh->Nelements*mesh->Nverts, 
	       MPI_DFLOAT, 0, 666, MPI_COMM_WORLD);
    }
    
    if(rank==0){
      MPI_Status status;
      MPI_Recv(tmpEX, allNelements[r]*mesh->Nverts, 
	       MPI_DFLOAT, r, 666, MPI_COMM_WORLD, &status);
      MPI_Recv(tmpEY, allNelements[r]*mesh->Nverts, 
	       MPI_DFLOAT, r, 666, MPI_COMM_WORLD, &status);
      
      for(iint e=0;e<allNelements[r];++e){
	fprintf(fp, "        ");
	for(iint n=0;n<mesh->Nverts;++n)
	  fprintf(fp, "%g %g 0\n", tmpEX[e*mesh->Nverts+n], tmpEY[e*mesh->Nverts+n]);
      }
    }
  }    
  if(rank==0){
    free(tmpEX);
    free(tmpEY);
    
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    
    // write out rank
    fprintf(fp, "      <CellData Scalars=\"scalars\">\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"element_rank\" Format=\"ascii\">\n");
    
    for(iint r=0;r<size;++r){
      for(iint e=0;e<allNelements[r];++e)
	fprintf(fp, "         %d\n", r);
    }
    
    fprintf(fp, "       </DataArray>\n");
    fprintf(fp, "     </CellData>\n");
    
    fprintf(fp, "    <Cells>\n");
    fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

    iint cnt=0;
    for(iint r=0;r<size;++r){
      for(iint e=0;e<allNelements[r];++e){
	for(iint n=0;n<mesh->Nverts;++n){
	  fprintf(fp, "%d ", cnt);
	  ++cnt;
	}
	fprintf(fp, "\n");
      }
    }

    fprintf(fp, "        </DataArray>\n");
    
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
    for(iint e=0;e<totalNelements;++e){
      if(e%10==0) fprintf(fp, "        ");
      fprintf(fp, "%d ", (e+1)*mesh->Nverts);
      if(((e+1)%10==0) || (e==totalNelements-1))
	fprintf(fp, "\n");
    }
    fprintf(fp, "       </DataArray>\n");
    
    fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
    for(iint e=0;e<totalNelements;++e){
      if(e%10==0) fprintf(fp, "        ");
      fprintf(fp, "5 ");
      if(((e+1)%10==0) || e==(totalNelements-1))
	fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
