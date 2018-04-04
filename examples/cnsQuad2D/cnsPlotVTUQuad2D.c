#include "cnsQuad2D.h"

// interpolate data to plot nodes and save to file (one per process
void cnsPlotVTUQuad2D(cns_t *cns, char *fileName){

  mesh2D *mesh = cns->mesh;
  
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
      dfloat plotxn = 0, plotyn = 0;

      for(int m=0;m<mesh->Np;++m){
        plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,0.);
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
        dfloat pm = mesh->q[e*mesh->Np*mesh->Nfields+m];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // write out velocity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"2\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0;
      for(int m=0;m<mesh->Np;++m){

	dfloat um = mesh->q[e*mesh->Np*mesh->Nfields+m+mesh->Np];
	dfloat vm = mesh->q[e*mesh->Np*mesh->Nfields+m+mesh->Np*2];
        //
        plotun += mesh->plotInterp[n*mesh->Np+m]*um;
        plotvn += mesh->plotInterp[n*mesh->Np+m]*vm;
      }
    
      fprintf(fp, "       ");
      fprintf(fp, "%g %g\n", plotun, plotvn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // write out vorticity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
  dfloat vort[mesh->plotNp];
  for(dlong e=0;e<mesh->Nelements;++e){

    for(int n=0;n<mesh->Np;++n){
      dfloat ur = 0, us = 0, vr = 0, vs = 0;

      for(int m=0;m<mesh->Np;++m){

	int qbase = e*mesh->Np*mesh->Nfields + m;

	dfloat rm = mesh->q[qbase+0*mesh->Np];
	dfloat um = mesh->q[qbase+1*mesh->Np]/rm;
	dfloat vm = mesh->q[qbase+2*mesh->Np]/rm;
	
	ur += mesh->Dr[n*mesh->Np+m]*um;
	us += mesh->Ds[n*mesh->Np+m]*um;
	vr += mesh->Dr[n*mesh->Np+m]*vm;
	vs += mesh->Ds[n*mesh->Np+m]*vm;
      }

      dfloat rx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + n + RXID*mesh->Np];
      dfloat sx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + n + SXID*mesh->Np];
      dfloat ry = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + n + RYID*mesh->Np];
      dfloat sy = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + n + SYID*mesh->Np];
      
      vort[n] =rx*vr + sx*vs - ry*ur - sy*us;
    }
    
    for(int n=0;n<mesh->plotNp;++n){

      dfloat plotvort = 0;
      for(int m=0;m<mesh->Np;++m){
        plotvort += mesh->plotInterp[n*mesh->Np+m]*vort[m];
      }
      
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotvort);
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
