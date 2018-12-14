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

#include "cds.h"

// interpolate data to plot nodes and save to file (one per process
void cdsPlotVTU(cds_t *cds, char *fileName){

  mesh_t *mesh = cds->mesh;
  
  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  FILE *fp;
  
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
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0,  plotzn = 0;

      for(int m=0;m<mesh->Np;++m){
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
  
  for(int fld=0; fld < cds->NSfields; fld++){
      // write out pressure
      fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Scalar%d\" Format=\"ascii\">\n",fld);
      for(dlong e=0;e<mesh->Nelements;++e){
	for(int n=0;n<mesh->plotNp;++n){
	  dfloat plotsn = 0.0;
	  for(int m=0;m<mesh->Np;++m){
	    dlong id = m+e*mesh->Np;
	    dfloat sm = cds->S[id+0*cds->sOffset];
	    plotsn += mesh->plotInterp[n*mesh->Np+m]*sm;
	  }

	  fprintf(fp, "       ");
	  fprintf(fp, "%g\n", plotsn);
	}
      }
      fprintf(fp, "       </DataArray>\n");
    }
  
	if (cds->dim==2) {
	  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"2\" Format=\"ascii\">\n");
	  for(dlong e=0;e<mesh->Nelements;++e){
	    for(int n=0;n<mesh->plotNp;++n){
	      dfloat plotun = 0, plotvn = 0;
	      for(int m=0;m<mesh->Np;++m){
		dlong id = m+e*mesh->Np;
		dfloat um = cds->U[id+0*offset];
		dfloat vm = cds->U[id+1*offset];

		plotun += mesh->plotInterp[n*mesh->Np+m]*um;
		plotvn += mesh->plotInterp[n*mesh->Np+m]*vm;
	      }
      
	      fprintf(fp, "       ");
	      fprintf(fp, "%g %g\n", plotun, plotvn);
	    }
	  }
	  fprintf(fp, "       </DataArray>\n");
	  fprintf(fp, "     </PointData>\n");
	} else {
	  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
	  for(dlong e=0;e<mesh->Nelements;++e){
	    for(int n=0;n<mesh->plotNp;++n){
	      dfloat plotun = 0, plotvn = 0, plotwn = 0;
	      for(int m=0;m<mesh->Np;++m){
		dlong id = m+e*mesh->Np;
		dfloat um = cds->U[id+0*offset];
		dfloat vm = cds->U[id+1*offset];
		dfloat wm = cds->U[id+2*offset];

		plotun += mesh->plotInterp[n*mesh->Np+m]*um;
		plotvn += mesh->plotInterp[n*mesh->Np+m]*vm;
		plotwn += mesh->plotInterp[n*mesh->Np+m]*wm;
	      }
      
	      fprintf(fp, "       ");
	      fprintf(fp, "%g %g %g\n", plotun, plotvn,plotwn);
	    }
	  }
	  fprintf(fp, "       </DataArray>\n");
	  fprintf(fp, "     </PointData>\n");
	}
  
	fprintf(fp, "    <Cells>\n");
      fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
      for(dlong e=0;e<mesh->Nelements;++e){
	for(int n=0;n<mesh->plotNelements;++n){
	  fprintf(fp, "       ");
	  for(int m=0;m<mesh->plotNverts;++m){
	    fprintf(fp, "%d ", e*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
	  }
	  fprintf(fp, "\n");
	}
      }
  
	fprintf(fp, "        </DataArray>\n");
  
      fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
      dlong cnt = 0;
      for(dlong e=0;e<mesh->Nelements;++e){
	for(int n=0;n<mesh->plotNelements;++n){
	  cnt += mesh->plotNverts;
	  fprintf(fp, "       ");
	  fprintf(fp, "%d\n", cnt);
	}
      }
	fprintf(fp, "       </DataArray>\n");
  
      fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh->Nelements;++e){
	for(int n=0;n<mesh->plotNelements;++n){
	  if(mesh->NfaceVertices==2)
	    fprintf(fp, "5\n");
	  else
	    fprintf(fp, "10\n");
	}
      }
	fprintf(fp, "        </DataArray>\n");
      fprintf(fp, "      </Cells>\n");
      fprintf(fp, "    </Piece>\n");
      fprintf(fp, "  </UnstructuredGrid>\n");
      fprintf(fp, "</VTKFile>\n");
      fclose(fp);

      }
