#include "bns.h"

// interpolate data to plot nodes and save to file (one per process
void bnsPlotVTU(bns_t *bns, char *fileName){

  mesh_t *mesh = bns->mesh;
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
      dfloat plotxn = 0, plotyn = 0, plotzn=0;

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
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(int m=0;m<mesh->Np;++m){
	        const dlong base = e*bns->Nfields*mesh->Np + m;
           dfloat rho = bns->q[base + 0*mesh->Np];
           dfloat pm  = bns->sqrtRT*bns->sqrtRT*rho; // need to be modified
          plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  


#if 0
  // calculate plot vorticity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"VorticityDivergence\" NumberOfComponents=\"2\" Format=\"ascii\">\n");
  dfloat *curlU = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *divU  = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat dUdr = 0, dUds = 0, dVdr = 0, dVds = 0;
      for(int m=0;m<mesh->Np;++m){
        int base = bns->Nfields*(m + e*mesh->Np);
        dfloat rho = bns->q[base + 0];
        dfloat u = bns->q[1 + base]*bns->sqrtRT/rho;
        dfloat v = bns->q[2 + base]*bns->sqrtRT/rho;
      	dUdr += mesh->Dr[n*mesh->Np+m]*u;
      	dUds += mesh->Ds[n*mesh->Np+m]*u;
      	dVdr += mesh->Dr[n*mesh->Np+m]*v;
      	dVds += mesh->Ds[n*mesh->Np+m]*v;
      }
      dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
      dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
      dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
      dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];

      dfloat dUdx = rx*dUdr + sx*dUds;
      dfloat dUdy = ry*dUdr + sy*dUds;
      dfloat dVdx = rx*dVdr + sx*dVds;
      dfloat dVdy = ry*dVdr + sy*dVds;
      
      curlU[n] = dVdx-dUdy;
      divU[n] = dUdx+dVdy;
    }
    
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotCurlUn = 0;
      dfloat plotDivUn = 0;
      for(int m=0;m<mesh->Np;++m){
        plotCurlUn += mesh->plotInterp[n*mesh->Np+m]*curlU[m];
        plotDivUn += mesh->plotInterp[n*mesh->Np+m]*divU[m];	
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g %g\n", plotCurlUn, plotDivUn);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  free(curlU);
  free(divU);
  #endif



  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0, plotwn=0;
      for(int m=0;m<mesh->Np;++m){
        dlong base = e*bns->Nfields*mesh->Np + m;
        dfloat rho = bns->q[base + 0*mesh->Np];
        dfloat um  = bns->q[base + 1*mesh->Np]*bns->sqrtRT/rho;
        dfloat vm  = bns->q[base + 2*mesh->Np]*bns->sqrtRT/rho;
        dfloat wm  = 0; 
        if(bns->dim==3)
          wm  = bns->q[base + 3*mesh->Np]*bns->sqrtRT/rho;
        //
        plotun += mesh->plotInterp[n*mesh->Np+m]*um;
        plotvn += mesh->plotInterp[n*mesh->Np+m]*vm;
        plotwn += mesh->plotInterp[n*mesh->Np+m]*wm;
        
      }
    
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotun, plotvn, plotwn);
    }
  }
  fprintf(fp, "       </DataArray>\n");




  
  fprintf(fp, "     </PointData>\n");
  
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
      if(bns->dim==2)
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
