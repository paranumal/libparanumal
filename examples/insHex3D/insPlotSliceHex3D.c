#include "insHex3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotSliceHex3D(ins_t *ins, char *fileName, const int Nslices, const char** dim, const dfloat* c){

  mesh3D *mesh = ins->mesh;
  
  //find number of sliced elements
  dfloat *coord;
  dlong NslicedElements = 0;
  dlong NslicedSubElements = 0;
  int *sliceFlag = (int *) calloc(mesh->Nelements,sizeof(int));
  int *sliceSubFlag = (int*) calloc(mesh->Nelements*mesh->plotNelements,sizeof(int));

  dfloat *plotCoord = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

  for (int i=0;i<Nslices;i++) {
    if (strstr(dim[i],"x")) coord = mesh->x;
    if (strstr(dim[i],"y")) coord = mesh->y;
    if (strstr(dim[i],"z")) coord = mesh->z;

    for (dlong e=0;e<mesh->Nelements;e++) {
      for(int n=0;n<mesh->plotNp;++n){
        plotCoord[n] = 0;
        for(int m=0;m<mesh->Np;++m){
          plotCoord[n] += mesh->plotInterp[n*mesh->Np+m]*coord[m+e*mesh->Np];
        }
      }

      for (int k=0;k<mesh->plotNelements;k++) {
        int id0 = mesh->plotEToV[k*mesh->plotNverts+0];
        int id1 = mesh->plotEToV[k*mesh->plotNverts+1];
        int id2 = mesh->plotEToV[k*mesh->plotNverts+2];
        int id3 = mesh->plotEToV[k*mesh->plotNverts+3];
        
        dfloat coord0 = plotCoord[id0];
        dfloat coord1 = plotCoord[id1];
        dfloat coord2 = plotCoord[id2];
        dfloat coord3 = plotCoord[id3];

        dfloat c0 = (c[i]-coord0)/(coord1-coord0);
        dfloat c1 = (c[i]-coord1)/(coord2-coord1);
        dfloat c2 = (c[i]-coord2)/(coord0-coord2);
        dfloat c3 = (c[i]-coord0)/(coord3-coord0);
        dfloat c4 = (c[i]-coord1)/(coord3-coord1);
        dfloat c5 = (c[i]-coord2)/(coord3-coord2);

        int cnt = 0;
        cnt += ((c0>=0.0)&&(c0<=1.0));
        cnt += ((c1>=0.0)&&(c1<=1.0));
        cnt += ((c2>=0.0)&&(c2<=1.0));
        cnt += ((c3>=0.0)&&(c3<=1.0));
        cnt += ((c4>=0.0)&&(c4<=1.0));
        cnt += ((c5>=0.0)&&(c5<=1.0));

        if (cnt) {
          if (sliceFlag[e]==0) NslicedElements++;
          sliceFlag[e] = 1;
          if (sliceSubFlag[e*mesh->plotNelements+k]==0) NslicedSubElements++;
          sliceSubFlag[e*mesh->plotNelements+k] = 1;
        }
      }
    }
  }

  //write sliced data to file
  FILE *fp;
  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", NslicedElements*mesh->plotNp, NslicedSubElements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
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
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotpn = 0;
      for(int m=0;m<mesh->Np;++m){
       const int offset = ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs); 
             const int id     = offset + m+e*mesh->Np;
        dfloat pm = ins->P[id];
        plotpn += mesh->plotInterp[n*mesh->Np+m]*pm;
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotpn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // write out divergence
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Divergence\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;

    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotDiv = 0;
      for(int m=0;m<mesh->Np;++m){
        int id = m+e*mesh->Np;
        plotDiv += mesh->plotInterp[n*mesh->Np+m]*ins->Div[id];
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotDiv);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // calculate plot vorticity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotVxn = 0, plotVyn = 0, plotVzn = 0 ;
      dfloat plotDivUn = 0;
      for(int m=0;m<mesh->Np;++m){
        int id = m+e*mesh->Np;
        plotVxn   += mesh->plotInterp[n*mesh->Np+m]*ins->Vx[id];
        plotVyn   += mesh->plotInterp[n*mesh->Np+m]*ins->Vy[id];
        plotVzn   += mesh->plotInterp[n*mesh->Np+m]*ins->Vz[id];
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotVxn, plotVyn, plotVzn);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotun = 0, plotvn = 0, plotwn=0;
      for(int m=0;m<mesh->Np;++m){
        int id = m+e*mesh->Np;
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
  dlong ccnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      if (sliceSubFlag[e*mesh->plotNelements+n]==0) continue;
      fprintf(fp, "       ");
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%d ", ccnt*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
    ccnt++;
  }
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  ccnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      if (sliceSubFlag[e*mesh->plotNelements+n]==0) continue;
      ccnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", ccnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(int e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(int n=0;n<mesh->plotNelements;++n){
      if (sliceSubFlag[e*mesh->plotNelements+n]==0) continue;
      fprintf(fp, "10\n"); // TET code ?
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  free(sliceFlag);
  free(sliceSubFlag);
}
