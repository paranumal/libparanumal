#include "insTet3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotContourTet3D(ins_t *ins, char *fileName, const char* options){

  mesh3D *mesh = ins->mesh;

  int Nlevels = 10;
  //dfloat levels[10] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
  dfloat levels[10] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
  dfloat tol = 1E-3;

  //if (strstr(options,"ADPATIVECONTOUR"))
    //meshPlotAdaptiveContour3D(mesh, fileName, Vort, Nlevels, levels, tol);
  //else 
    //meshPlotContour3D(mesh, fileName, Vort, Nlevels, levels);

  int *plotFlag = (int*) calloc(mesh->Nelements,sizeof(int));
  int *plotSubFlag = (int*) calloc(mesh->Nelements*mesh->plotNelements,sizeof(int));
  dfloat *plotvx = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
  dfloat *plotvy = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
  dfloat *plotvz = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

  dlong NcontourElements =0;
  dlong plotElements =0;

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNp;++n){
      plotvx[n] = 0; plotvy[n] = 0; plotvz[n] = 0;
      for(int m=0;m<mesh->Np;++m){
        plotvx[n] += mesh->plotInterp[n*mesh->Np+m]*ins->Vx[m+e*mesh->Np];
        plotvy[n] += mesh->plotInterp[n*mesh->Np+m]*ins->Vy[m+e*mesh->Np];
        plotvz[n] += mesh->plotInterp[n*mesh->Np+m]*ins->Vz[m+e*mesh->Np];
      }
    }

    for (int k=0;k<mesh->plotNelements;k++) {
      int id0 = mesh->plotEToV[k*mesh->plotNverts+0];
      int id1 = mesh->plotEToV[k*mesh->plotNverts+1];
      int id2 = mesh->plotEToV[k*mesh->plotNverts+2];
      int id3 = mesh->plotEToV[k*mesh->plotNverts+3];

      dfloat umin = sqrt(plotvx[id0]*plotvx[id0] + plotvy[id0]*plotvy[id0] + plotvz[id0]*plotvz[id0]);
      dfloat umax = sqrt(plotvx[id0]*plotvx[id0] + plotvy[id0]*plotvy[id0] + plotvz[id0]*plotvz[id0]);  
      umin = mymin(umin, sqrt(plotvx[id1]*plotvx[id1] + plotvy[id1]*plotvy[id1] + plotvz[id1]*plotvz[id1]));
      umax = mymax(umax, sqrt(plotvx[id1]*plotvx[id1] + plotvy[id1]*plotvy[id1] + plotvz[id1]*plotvz[id1]));
      umin = mymin(umin, sqrt(plotvx[id2]*plotvx[id2] + plotvy[id2]*plotvy[id2] + plotvz[id2]*plotvz[id2]));
      umax = mymax(umax, sqrt(plotvx[id2]*plotvx[id2] + plotvy[id2]*plotvy[id2] + plotvz[id2]*plotvz[id2]));
      umin = mymin(umin, sqrt(plotvx[id3]*plotvx[id3] + plotvy[id3]*plotvy[id3] + plotvz[id3]*plotvz[id3]));
      umax = mymax(umax, sqrt(plotvx[id3]*plotvx[id3] + plotvy[id3]*plotvy[id3] + plotvz[id3]*plotvz[id3]));

      for (int lev=0;lev<Nlevels;lev++){
        if((umin<=levels[lev]) && (umax>=levels[lev])){
          NcontourElements++;
          if (plotFlag[e]==0) plotElements++;
          plotFlag[e] = 1;
          plotSubFlag[e*mesh->plotNelements+k] = 1;
          break;
        }  
      }
    }
  }
  free(plotvx); free(plotvy); free(plotvz);

  FILE *fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
          plotElements*mesh->plotNp, 
          NcontourElements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
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
  
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotvx = 0, plotvy = 0, plotvz = 0;
      for(int m=0;m<mesh->Np;++m){
        plotvx += mesh->plotInterp[n*mesh->Np+m]*ins->Vx[m+e*mesh->Np];
        plotvy += mesh->plotInterp[n*mesh->Np+m]*ins->Vy[m+e*mesh->Np];
        plotvz += mesh->plotInterp[n*mesh->Np+m]*ins->Vz[m+e*mesh->Np];
      }
      fprintf(fp, "       "); fprintf(fp, "%g\n", sqrt(plotvx*plotvx+plotvy*plotvy+plotvz*plotvz));
    }
  }
  fprintf(fp, "       </DataArray>\n");
  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  dlong cnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%d ", cnt*mesh->plotNp + mesh->plotEToV[k*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
    cnt++;
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  cnt=0;
  for(dlong e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      cnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (plotFlag[e]==0) continue;
    for(int k=0;k<mesh->plotNelements;++k){
      if (plotSubFlag[e*mesh->plotNelements+k]==0) continue;
      fprintf(fp, "10\n"); // TET code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  free(plotFlag);
  free(plotSubFlag);
}
