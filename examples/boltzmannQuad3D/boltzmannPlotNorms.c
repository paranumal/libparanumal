#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "boltzmannQuad3D.h"

// interpolate data to plot nodes and save to file (one per process
void boltzmannPlotNorms(mesh_t *mesh, char *fileNameBase, int tstep,dfloat *q){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint plotnp = (mesh->Np + mesh->Nq)/2;
  
  FILE *fp;
  char fileName[BUFSIZ];
  sprintf(fileName, "%s_%04d_%04d.vtu", fileNameBase, rank, tstep);
  
  printf("FILE = %s\n", fileName);

  fp = fopen(fileName, "w");
  
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	  mesh->Nelements*2*plotnp,
	  mesh->Nelements*2);

  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=0;j<mesh->Nq-i;++j){
	dfloat plotxn = mesh->x[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotyn = mesh->y[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotzn = mesh->z[j+i*mesh->Nq+e*mesh->Np];
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
      }
    }
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=mesh->Nq-i-1;j<mesh->Nq;++j){
	dfloat plotxn = mesh->x[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotyn = mesh->y[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotzn = mesh->z[j+i*mesh->Nq+e*mesh->Np];
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"mrab_levels\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np+mesh->Nq;++n){
      dfloat plotpn = mesh->MRABlevel[e];
      //      plotpn = plotpn*mesh->rho; // Get Pressure
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }

  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"q_value\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=0; j<mesh->Nq-i;++j){
	dfloat plotpn = q[j+i*mesh->Nq+1*mesh->Np+e*mesh->Np*mesh->Nfields];
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotpn);
      }
    }
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=mesh->Nq-i-1; j<mesh->Nq;++j){
	dfloat plotpn = q[j+i*mesh->Nq+1*mesh->Np+e*mesh->Np*mesh->Nfields];
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotpn);
      }
    }
  }
  fprintf(fp, "       </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint i=0;i<mesh->Nq;++i) {
      for(iint j=0;j<mesh->Nq-i;++j) {
	int id = j+i*mesh->Nq + 0*mesh->Np + mesh->Nfields*mesh->Np*e;
        dfloat plotpn = mesh->q[id];
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotpn);
      }
    }
    for(iint i=0;i<mesh->Nq;++i) {
      for(iint j=mesh->Nq-i-1;j<mesh->Nq;++j) {
	int id = j+i*mesh->Nq + 0*mesh->Np + mesh->Nfields*mesh->Np*e;
        dfloat plotpn = mesh->q[id];
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotpn);
      }
    }
  }
  fprintf(fp, "       </DataArray>\n");


  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=0;j<mesh->Nq-i;++j){
	int uid = j+i*mesh->Nq + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = j+i*mesh->Nq + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = j+i*mesh->Nq + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
	
	dfloat plotun = mesh->q[uid];
	dfloat plotvn = mesh->q[vid];
	dfloat plotwn = mesh->q[vid];
	
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotun, plotvn, plotwn);
      }
    }
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=mesh->Nq-i-1;j<mesh->Nq;++j){
	int uid = j+i*mesh->Nq + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = j+i*mesh->Nq + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = j+i*mesh->Nq + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
	
	dfloat plotun = mesh->q[uid];
	dfloat plotvn = mesh->q[vid];
	dfloat plotwn = mesh->q[vid];
	
	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotun, plotvn, plotwn);
      }
    }
  }
  fprintf(fp, "       </DataArray>\n");

  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
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
          
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=0;j<mesh->Nq-i;++j){
	int uid = j+i*mesh->Nq + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = j+i*mesh->Nq + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = j+i*mesh->Nq + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
	
	dfloat plotvort1n = vort[j+i*mesh->Nq+0*mesh->Np];
	dfloat plotvort2n = vort[j+i*mesh->Nq+1*mesh->Np];
	dfloat plotvort3n = vort[j+i*mesh->Nq+2*mesh->Np];

	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotvort1n, plotvort2n, plotvort3n);

	dfloat plotxn = mesh->x[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotyn = mesh->y[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotzn = mesh->z[j+i*mesh->Nq+e*mesh->Np];
	
	plotVortRadial[j+i*mesh->Nq + mesh->plotNp*e] =
	  plotxn*plotvort1n + 
	  plotyn*plotvort2n + 
	  plotzn*plotvort3n;
      }
    }

    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=mesh->Nq-i-1;j<mesh->Nq;++j){
	int uid = j+i*mesh->Nq + 1*mesh->Np + mesh->Nfields*mesh->Np*e;
	int vid = j+i*mesh->Nq + 2*mesh->Np + mesh->Nfields*mesh->Np*e;
	int wid = j+i*mesh->Nq + 3*mesh->Np + mesh->Nfields*mesh->Np*e;
	
	dfloat plotvort1n = vort[j+i*mesh->Nq+0*mesh->Np];
	dfloat plotvort2n = vort[j+i*mesh->Nq+1*mesh->Np];
	dfloat plotvort3n = vort[j+i*mesh->Nq+2*mesh->Np];

	fprintf(fp, "       ");
	fprintf(fp, "%g %g %g\n", plotvort1n, plotvort2n, plotvort3n);

	dfloat plotxn = mesh->x[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotyn = mesh->y[j+i*mesh->Nq+e*mesh->Np];
	dfloat plotzn = mesh->z[j+i*mesh->Nq+e*mesh->Np];
	
	plotVortRadial[j+i*mesh->Nq + mesh->plotNp*e] =
	  plotxn*plotvort1n + 
	  plotyn*plotvort2n + 
	  plotzn*plotvort3n;
      }
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"RadialVorticity\" NumberOfComponents=\"1\" Format=\"ascii\">\n");

  for(iint e=0;e<mesh->Nelements;++e){

    // compute vorticity
    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=0;j<mesh->Nq-i;++j){
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotVortRadial[j+i*mesh->Nq+e*mesh->plotNp]);
      }
    }

    for(iint i=0;i<mesh->Nq;++i){
      for(iint j=mesh->Nq-i-1;j<mesh->Nq;++j){
	fprintf(fp, "       ");
	fprintf(fp, "%g\n", plotVortRadial[j+i*mesh->Nq+e*mesh->plotNp]);
      }
    }
  }

  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "    </PointData>\n");
  
  free(vort);
  free(plotVortRadial);
    
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<mesh->Nelements;++e){
    fprintf(fp, "       ");
    fprintf(fp, "%d %d %d \n", 2*e*(mesh->Np + mesh->Nq)/2 + 0,2*e*(mesh->Np + mesh->Nq)/2 + mesh->Nq - 1,(2*e+1)*(mesh->Np+mesh->Nq)/2-1);
    fprintf(fp, "%d %d %d \n", (2*e+1)*(mesh->Np + mesh->Nq)/2 + 0,(2*e+2)*(mesh->Np + mesh->Nq)/2 - mesh->Nq,(2*e+2)*(mesh->Np+mesh->Nq)/2-1);
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<2*mesh->Nelements;++e){
      cnt += 3;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<2*mesh->Nelements;++e){
      fprintf(fp, "5\n");
  }

  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}
