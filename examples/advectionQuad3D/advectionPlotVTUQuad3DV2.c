#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "advectionQuad3D.h"

// interpolate data to plot nodes and save to file (one per process
void advectionPlotVTUQuad3DV2(mesh_t *mesh, char *fileNameBase, int tstep){

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
  
  
  // write out pressure
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");



  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(iint m=0;m<mesh->Np;++m){
	int id = m + 0*mesh->Np + mesh->Nfields*mesh->Np*e;
        dfloat pm = mesh->q[id];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }
      //
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");


  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0, plotwn = 0;
      for(iint m=0;m<mesh->Np;++m){
	int uid = m + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = m + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = m + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
		
        dfloat um = mesh->q[uid];
        dfloat vm = mesh->q[vid];
	dfloat wm = mesh->q[wid];
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

  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticit\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  dfloat *vort = (dfloat*) calloc(mesh->Np*mesh->dim, sizeof(dfloat));
  dfloat *plotVortRadial = (dfloat*) calloc(mesh->plotNp*mesh->Nelements, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

    // compute vorticity
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){

	int n = i+mesh->Nq*j;
	iint gbase = e*mesh->Np*mesh->Nvgeo + n;
	dfloat rx = mesh->vgeo[gbase+mesh->Np*RXID];
	dfloat sx = mesh->vgeo[gbase+mesh->Np*SXID];
	dfloat tx = mesh->vgeo[gbase+mesh->Np*TXID];
	
	dfloat ry = mesh->vgeo[gbase+mesh->Np*RYID];
	dfloat sy = mesh->vgeo[gbase+mesh->Np*SYID];
	dfloat ty = mesh->vgeo[gbase+mesh->Np*TYID];
	
	dfloat rz = mesh->vgeo[gbase+mesh->Np*RZID];
	dfloat sz = mesh->vgeo[gbase+mesh->Np*SZID];
	dfloat tz = mesh->vgeo[gbase+mesh->Np*TZID];
	
	dfloat dudr = 0, duds = 0;
	dfloat dvdr = 0, dvds = 0;
	dfloat dwdr = 0, dwds = 0;
	
	for(iint m=0;m<mesh->Nq;++m){
	  int basejm = m + j*mesh->Nq + mesh->Nfields*mesh->Np*e;
	  int basemi = i + m*mesh->Nq + mesh->Nfields*mesh->Np*e;
	  dfloat Dim = mesh->D[i*mesh->Nq+m];
	  dfloat Djm = mesh->D[j*mesh->Nq+m];
	  dudr += Dim*mesh->q[basejm + 1*mesh->Np]/mesh->q[basejm];
	  duds += Djm*mesh->q[basemi + 1*mesh->Np]/mesh->q[basemi];
	  dvdr += Dim*mesh->q[basejm + 2*mesh->Np]/mesh->q[basejm];
	  dvds += Djm*mesh->q[basemi + 2*mesh->Np]/mesh->q[basemi];
	  dwdr += Dim*mesh->q[basejm + 3*mesh->Np]/mesh->q[basejm];
	  dwds += Djm*mesh->q[basemi + 3*mesh->Np]/mesh->q[basemi];
	}
	int base = i + j*mesh->Nq + e*mesh->Np*mesh->Nfields;
	dfloat dudx = rx*dudr + sx*duds + tx*mesh->q[base + 1*mesh->Np]/mesh->q[base];
	dfloat dudy = ry*dudr + sy*duds + ty*mesh->q[base + 1*mesh->Np]/mesh->q[base];
	dfloat dudz = rz*dudr + sz*duds + tz*mesh->q[base + 1*mesh->Np]/mesh->q[base];

	dfloat dvdx = rx*dvdr + sx*dvds + tx*mesh->q[base + 2*mesh->Np]/mesh->q[base];
	dfloat dvdy = ry*dvdr + sy*dvds + ty*mesh->q[base + 2*mesh->Np]/mesh->q[base];
	dfloat dvdz = rz*dvdr + sz*dvds + tz*mesh->q[base + 2*mesh->Np]/mesh->q[base];

	dfloat dwdx = rx*dwdr + sx*dwds + tx*mesh->q[base + 3*mesh->Np]/mesh->q[base];
	dfloat dwdy = ry*dwdr + sy*dwds + ty*mesh->q[base + 3*mesh->Np]/mesh->q[base];
	dfloat dwdz = rz*dwdr + sz*dwds + tz*mesh->q[base + 3*mesh->Np]/mesh->q[base];

	base = i + j*mesh->Nq;
	vort[base+0*mesh->Np] = dwdy - dvdz;
	vort[base+1*mesh->Np] = dudz - dwdx;
	vort[base+2*mesh->Np] = dvdx - dudy;

      }
    }
    
    for(iint n=0;n<mesh->plotNp;++n){
      
      dfloat plotvort1n = 0, plotvort2n = 0, plotvort3n = 0;
      for(iint m=0;m<mesh->Np;++m){
	int uid = m + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = m + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = m + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
		
        dfloat vort1m = vort[m+0*mesh->Np];
	dfloat vort2m = vort[m+1*mesh->Np];
	dfloat vort3m = vort[m+2*mesh->Np];
        //
        plotvort1n += mesh->plotInterp[n*mesh->Np+m]*vort1m;
        plotvort2n += mesh->plotInterp[n*mesh->Np+m]*vort2m;
	plotvort3n += mesh->plotInterp[n*mesh->Np+m]*vort3m;
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotvort1n, plotvort2n, plotvort3n);

      dfloat plotxn = 0, plotyn = 0, plotzn = 0;
      for(iint m=0;m<mesh->Np;++m){
	plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
	plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
	plotzn += mesh->plotInterp[n*mesh->Np+m]*mesh->z[m+e*mesh->Np];
      }
      
      plotVortRadial[n + mesh->plotNp*e] =
	plotxn*plotvort1n + 
	plotyn*plotvort2n + 
	plotzn*plotvort3n;
      
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"RadialVorticity\" NumberOfComponents=\"1\" Format=\"ascii\">\n");

  for(iint e=0;e<mesh->Nelements;++e){

    // compute vorticity
    for(iint n=0;n<mesh->plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotVortRadial[n+e*mesh->plotNp]);
    }
  }

  fprintf(fp, "       </DataArray>\n");
  
  free(vort);
  free(plotVortRadial);

  // fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  // for(iint e=0;e<mesh->Nelements;++e){
  //   for(iint n=0;n<mesh->plotNp;++n){
  //     dfloat plotwxn = 0, plotwyn = 0, plotvn = 0, plotwzn = 0;
  //     for(iint m=0;m<mesh->Np;++m){
  //       dfloat wx = mesh->q[4 + mesh->Nfields*(m+e*mesh->Np)];
  //       dfloat wy = mesh->q[5 + mesh->Nfields*(m+e*mesh->Np)];
  //       dfloat wz = mesh->q[6 + mesh->Nfields*(m+e*mesh->Np)];
  //       //
  //       plotwxn += mesh->plotInterp[n*mesh->Np+m]*wx;
  //       plotwyn += mesh->plotInterp[n*mesh->Np+m]*wy;
  //       plotwzn += mesh->plotInterp[n*mesh->Np+m]*wz;
  //     }
      
  //     fprintf(fp, "       ");
  //     fprintf(fp, "%g %g %g\n", plotwxn, plotwyn,plotwzn);
  //   }
  // }
  // fprintf(fp, "       </DataArray>\n");



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
