#include "gradient.h"

// interpolate data to plot nodes and save to file (one per process
void gradientPlotVTU(gradient_t *gradient, char *fileName){

  mesh_t *mesh = gradient->mesh;

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
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;

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
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(int m=0;m<mesh->Np;++m){
        dfloat pm = mesh->q[e*mesh->Np*mesh->Nfields+m];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // write out velocity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0, plotwn = 0;
      for(int m=0;m<mesh->Np;++m){
        dfloat rm = mesh->q[e*mesh->Np*mesh->Nfields+m           ];
        dfloat um = mesh->q[e*mesh->Np*mesh->Nfields+m+mesh->Np  ]/rm;
        dfloat vm = mesh->q[e*mesh->Np*mesh->Nfields+m+mesh->Np*2]/rm;
        dfloat wm = mesh->q[e*mesh->Np*mesh->Nfields+m+mesh->Np*3]/rm;
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

  // write out vorticity (need to fix for 3D vorticity)
  const int Ncomponents = (gradient->dim==2) ? 1: 3;
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Gradient\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ncomponents);

  if(gradient->dim==3){
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->plotNp;++n){
	dfloat plotgradientqx = 0, plotgradientqy = 0, plotgradientqz = 0;
	for(int m=0;m<mesh->Np;++m){
	  dlong id = m+e*mesh->Np*3;
	  dfloat vortx = gradient->gradientq[id];
	  dfloat vorty = gradient->gradientq[id+mesh->Np];
	  dfloat vortz = gradient->gradientq[id+2*mesh->Np];
	  plotgradientqx += mesh->plotInterp[n*mesh->Np+m]*vortx;
        plotgradientqy += mesh->plotInterp[n*mesh->Np+m]*vorty;
        plotgradientqz += mesh->plotInterp[n*mesh->Np+m]*vortz;
	}
	
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotgradientqx, plotgradientqy, plotgradientqz);
      }
    }
  }
  else{
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->plotNp;++n){
	dfloat plotgradientq = 0;
	for(int m=0;m<mesh->Np;++m){
	  dlong id = m+e*mesh->Np;
	  dfloat vort = gradient->gradientq[id];
	  plotgradientq += mesh->plotInterp[n*mesh->Np+m]*vort;
	}
	
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotgradientq);
      }
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
      if(gradient->dim==2)
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
