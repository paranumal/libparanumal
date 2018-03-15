#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "ins3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotVTU3D(ins_t *ins, char *fileName){

  mesh3D *mesh = ins->mesh;
  
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


  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(iint m=0;m<mesh->Np;++m){
       const iint offset = ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs); 
	     const iint id     = offset + m+e*mesh->Np;
        dfloat pm = ins->P[id];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  // fprintf(fp, "       </DataArray>\n");

  // // calculate plot vorticity
  // fprintf(fp, "        <DataArray type=\"Float32\" Name=\"VorticityDivergence\" NumberOfComponents=\"4\" Format=\"ascii\">\n");
  // dfloat *Vx   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  // dfloat *Vy   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  // dfloat *Vz   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  // dfloat *divU = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  
  // for(iint e=0;e<mesh->Nelements;++e){
  //   for(iint n=0;n<mesh->Np;++n){
  //     dfloat dUdr = 0, dUds = 0, dUdt = 0 ;
  //     dfloat dVdr = 0, dVds = 0, dVdt = 0 ;
  //     dfloat dWdr = 0, dWds = 0, dWdt = 0 ; 
  //     for(iint m=0;m<mesh->Np;++m){
  //     	iint id = m+e*mesh->Np;
  //     	id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

  //     	dUdr += mesh->Dr[n*mesh->Np+m]*ins->U[id];
  //       dUds += mesh->Ds[n*mesh->Np+m]*ins->U[id];
  //       dUdt += mesh->Dt[n*mesh->Np+m]*ins->U[id];

  //     	dVdr += mesh->Dr[n*mesh->Np+m]*ins->V[id];
  //       dVds += mesh->Ds[n*mesh->Np+m]*ins->V[id];
  //       dVdt += mesh->Dt[n*mesh->Np+m]*ins->V[id];

  //       dWdr += mesh->Dr[n*mesh->Np+m]*ins->W[id];
  //       dWds += mesh->Ds[n*mesh->Np+m]*ins->W[id];
  //       dWdt += mesh->Dt[n*mesh->Np+m]*ins->W[id];
  //     }

  //     dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
  //     dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
  //     dfloat rz = mesh->vgeo[e*mesh->Nvgeo+RZID];    
      
  //     dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
  //     dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];
  //     dfloat sz = mesh->vgeo[e*mesh->Nvgeo+SZID];    
     
  //     dfloat tx = mesh->vgeo[e*mesh->Nvgeo+TXID];
  //     dfloat ty = mesh->vgeo[e*mesh->Nvgeo+TYID];
  //     dfloat tz = mesh->vgeo[e*mesh->Nvgeo+TZID];    

  //     dfloat dUdx = rx*dUdr + sx*dUds + tx*dUdt;
  //     dfloat dUdy = ry*dUdr + sy*dUds + ty*dUdt;
  //     dfloat dUdz = rz*dUdr + sz*dUds + tz*dUdt;
    
  //     dfloat dVdx = rx*dVdr + sx*dVds + tx*dVdt;
  //     dfloat dVdy = ry*dVdr + sy*dVds + ty*dVdt;
  //     dfloat dVdz = rz*dVdr + sz*dVds + tz*dVdt;
      
  //     dfloat dWdx = rx*dWdr + sx*dWds + tx*dWdt;
  //     dfloat dWdy = ry*dWdr + sy*dWds + ty*dWdt;
  //     dfloat dWdz = rz*dWdr + sz*dWds + tz*dWdt;
      
  //     // Compute vorticity Vector
  //     Vx[n] = dWdy-dVdz;
  //     Vy[n] = dUdz-dWdx;
  //     Vz[n] = dVdx-dUdy;
  //     //
  //     divU[n] = dUdx + dVdy + dWdz; 
  //   }
    
  //   for(iint n=0;n<mesh->plotNp;++n){
  //     dfloat plotVxn = 0, plotVyn = 0, plotVzn = 0 ;
  //     dfloat plotDivUn = 0;
  //     for(iint m=0;m<mesh->Np;++m){
  //       plotVxn   += mesh->plotInterp[n*mesh->Np+m]*Vx[m];
  //       plotVyn   += mesh->plotInterp[n*mesh->Np+m]*Vy[m];
  //       plotVzn   += mesh->plotInterp[n*mesh->Np+m]*Vz[m];
  //       plotDivUn += mesh->plotInterp[n*mesh->Np+m]*divU[m];	
  //     }
  //     fprintf(fp, "       ");
  //     fprintf(fp, "%g %g %g %g\n", plotVxn, plotVyn, plotVzn, plotDivUn);
  //   }
  // }
  // free(Vx);
  // free(Vy);
  // free(Vz);
  // free(divU);

  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0, plotwn=0;
      for(iint m=0;m<mesh->Np;++m){
	iint id = m+e*mesh->Np;
	id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
        dfloat um = ins->U[id];
        dfloat vm = ins->V[id];
        dfloat wm = ins->W[id];
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
