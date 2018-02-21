#include "ins3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotSlice3D(ins_t *ins, char *fileName, const int Nslices, const char** dim, const dfloat* c){

  mesh3D *mesh = ins->mesh;
  
  //find number of sliced elements
  dfloat *coord;
  iint NslicedElements = 0;
  iint NslicedSubElements = 0;
  int *sliceFlag = (int *) calloc(mesh->Nelements,sizeof(int));
  int *sliceSubFlag = (int*) calloc(mesh->Nelements*mesh->plotNelements,sizeof(int));

  dfloat *plotCoord = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

  for (int i=0;i<Nslices;i++) {
    if (strstr(dim[i],"x")) coord = mesh->x;
    if (strstr(dim[i],"y")) coord = mesh->y;
    if (strstr(dim[i],"z")) coord = mesh->z;

    for (iint e=0;e<mesh->Nelements;e++) {
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


#if 1
  //write sliced data to file
  FILE *fp;
  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", NslicedElements*mesh->plotNp, NslicedSubElements);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
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
    if (sliceFlag[e]==0) continue;
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
  fprintf(fp, "       </DataArray>\n");

  // write out divergence
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Divergence\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;

    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotDiv = 0;
      for(iint m=0;m<mesh->Np;++m){
        iint id = m+e*mesh->Np;
        plotDiv += mesh->plotInterp[n*mesh->Np+m]*ins->Div[id];
      }
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotDiv);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  // calculate plot vorticity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotVxn = 0, plotVyn = 0, plotVzn = 0 ;
      dfloat plotDivUn = 0;
      for(iint m=0;m<mesh->Np;++m){
        iint id = m+e*mesh->Np;
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
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
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
  iint ccnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(iint n=0;n<mesh->plotNelements;++n){
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
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(iint n=0;n<mesh->plotNelements;++n){
      if (sliceSubFlag[e*mesh->plotNelements+n]==0) continue;
      ccnt += mesh->plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", ccnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<mesh->Nelements;++e){
    if (sliceFlag[e]==0) continue;
    for(iint n=0;n<mesh->plotNelements;++n){
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
#else
    int *elementType = (int *) calloc(NslicedElements*mesh->plotNelements,sizeof(int));

    dfloat *plotx = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *ploty = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotz = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *xn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *yn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *zn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

    dfloat *plotU = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotV = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotW = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotP = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *Un = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Vn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Wn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Pn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

    dfloat *plotVortx = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotVorty = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotVortz = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));
    dfloat *plotDiv = (dfloat *) calloc(NslicedElements*mesh->plotNelements*3,sizeof(dfloat));

    dfloat *Vortx = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
    dfloat *Vorty = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
    dfloat *Vortz = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
    dfloat *Div = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
    dfloat *Vortxn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Vortyn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Vortzn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));
    dfloat *Divn = (dfloat *) calloc(mesh->plotNp,sizeof(dfloat));

    iint cnt = 0;
    NslicedElements =0;
    for (iint e=0;e<mesh->Nelements;e++) {
      if (sliceFlag[e]==0) continue;

      //compute div and vort
      for(iint n=0;n<mesh->Np;++n){
        dfloat dUdr = 0, dUds = 0, dUdt = 0 ;
        dfloat dVdr = 0, dVds = 0, dVdt = 0 ;
        dfloat dWdr = 0, dWds = 0, dWdt = 0 ; 
        for(iint m=0;m<mesh->Np;++m){
          iint id = m+e*mesh->Np;
          id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

          dUdr += mesh->Dr[n*mesh->Np+m]*ins->U[id];
          dUds += mesh->Ds[n*mesh->Np+m]*ins->U[id];
          dUdt += mesh->Dt[n*mesh->Np+m]*ins->U[id];

          dVdr += mesh->Dr[n*mesh->Np+m]*ins->V[id];
          dVds += mesh->Ds[n*mesh->Np+m]*ins->V[id];
          dVdt += mesh->Dt[n*mesh->Np+m]*ins->V[id];

          dWdr += mesh->Dr[n*mesh->Np+m]*ins->W[id];
          dWds += mesh->Ds[n*mesh->Np+m]*ins->W[id];
          dWdt += mesh->Dt[n*mesh->Np+m]*ins->W[id];
        }

        dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
        dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
        dfloat rz = mesh->vgeo[e*mesh->Nvgeo+RZID];    
        
        dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
        dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];
        dfloat sz = mesh->vgeo[e*mesh->Nvgeo+SZID];    
       
        dfloat tx = mesh->vgeo[e*mesh->Nvgeo+TXID];
        dfloat ty = mesh->vgeo[e*mesh->Nvgeo+TYID];
        dfloat tz = mesh->vgeo[e*mesh->Nvgeo+TZID];    

        dfloat dUdx = rx*dUdr + sx*dUds + tx*dUdt;
        dfloat dUdy = ry*dUdr + sy*dUds + ty*dUdt;
        dfloat dUdz = rz*dUdr + sz*dUds + tz*dUdt;
      
        dfloat dVdx = rx*dVdr + sx*dVds + tx*dVdt;
        dfloat dVdy = ry*dVdr + sy*dVds + ty*dVdt;
        dfloat dVdz = rz*dVdr + sz*dVds + tz*dVdt;
        
        dfloat dWdx = rx*dWdr + sx*dWds + tx*dWdt;
        dfloat dWdy = ry*dWdr + sy*dWds + ty*dWdt;
        dfloat dWdz = rz*dWdr + sz*dWds + tz*dWdt;
        
        // Compute vorticity Vector
        Vortx[n] = dWdy-dVdz;
        Vorty[n] = dUdz-dWdx;
        Vortz[n] = dVdx-dUdy;
        
        Div[n] = dUdx + dVdy + dWdz; 
      }
  
      //interpolate 3D element to plot nodes
      for(int n=0;n<mesh->plotNp;++n){  
        xn[n] = 0; yn[n] = 0; zn[n] = 0;
        Un[n] = 0.; Vn[n] = 0.; Wn[n] = 0.; Pn[n] = 0.;
        Vortxn[n] = 0.; Vortyn[n] = 0.; Vortzn[n] = 0.; Divn[n] = 0.;

        for(int m=0;m<mesh->Np;++m){
          iint id = m+e*mesh->Np;
          id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

          xn[n] += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
          yn[n] += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
          zn[n] += mesh->plotInterp[n*mesh->Np+m]*mesh->z[m+e*mesh->Np];
          Un[n] += mesh->plotInterp[n*mesh->Np+m]*ins->U[id];
          Vn[n] += mesh->plotInterp[n*mesh->Np+m]*ins->V[id];
          Wn[n] += mesh->plotInterp[n*mesh->Np+m]*ins->W[id];
          Pn[n] += mesh->plotInterp[n*mesh->Np+m]*ins->P[id];
          Vortxn[n] += mesh->plotInterp[n*mesh->Np+m]*Vortx[m];
          Vortyn[n] += mesh->plotInterp[n*mesh->Np+m]*Vorty[m];
          Vortzn[n] += mesh->plotInterp[n*mesh->Np+m]*Vortz[m];
          Divn[n] += mesh->plotInterp[n*mesh->Np+m]*Div[m];
        }
      }

      //check if a submesh element is sliced
      for (int n=0;n<mesh->plotNelements;n++) {
        int id0 = mesh->plotEToV[n*mesh->Nverts+0];
        int id1 = mesh->plotEToV[n*mesh->Nverts+1];
        int id2 = mesh->plotEToV[n*mesh->Nverts+2];
        int id3 = mesh->plotEToV[n*mesh->Nverts+3];
              
        //edge vertices
        int cid0[6], cid1[6];
        cid0[0] = id0; cid1[0] = id1;
        cid0[1] = id1; cid1[1] = id2;
        cid0[2] = id2; cid1[2] = id0;
        cid0[3] = id2; cid1[3] = id3;
        cid0[4] = id1; cid1[4] = id3;
        cid0[5] = id0; cid1[5] = id3;

        for (int k=0;k<Nslices;k++) {
          if (strstr(dim[k],"x")) coord = xn;
          if (strstr(dim[k],"y")) coord = yn;
          if (strstr(dim[k],"z")) coord = zn;

          dfloat cc[6];
          for (int i=0;i<6;i++) {
            int id0 = cid0[i], id1 = cid1[i];
            cc[i] = (c[k]-coord[id0])/(coord[id1]-coord[id0]);
          }
          int ccnt[7];
          ccnt[0] = ((cc[0]>=0.0)&&(cc[0]<1.0));
          ccnt[1] = ((cc[1]>=0.0)&&(cc[1]<1.0));
          ccnt[2] = ((cc[2]>=0.0)&&(cc[2]<1.0));
          ccnt[3] = ((cc[3]>0.0) &&(cc[3]<=1.0));
          ccnt[4] = ((cc[4]>0.0) &&(cc[4]<1.0));
          ccnt[5] = ((cc[5]>0.0) &&(cc[5]<1.0));

          ccnt[6] = 0;
          for (int i=0;i<6;i++) ccnt[6] += ccnt[i];

          if (ccnt[6]>4) printf("huh, again?\n");
          if ((ccnt[6]==3)){ //if sliced, add the data to the array
            elementType[NslicedElements++] = ccnt[6];
            for (int i=0;i<6;i++) {
              if (ccnt[i]) {
                int id0 = cid0[i], id1 = cid1[i];
                plotx[cnt] = (1-cc[i])*xn[id0] + cc[i]*xn[id1];
                ploty[cnt] = (1-cc[i])*yn[id0] + cc[i]*yn[id1];
                plotz[cnt] = (1-cc[i])*zn[id0] + cc[i]*zn[id1];
                
                plotU[cnt] = (1-cc[i])*Un[id0] + cc[i]*Un[id1];
                plotV[cnt] = (1-cc[i])*Vn[id0] + cc[i]*Vn[id1];
                plotW[cnt] = (1-cc[i])*Wn[id0] + cc[i]*Wn[id1];
                plotP[cnt] = (1-cc[i])*Pn[id0] + cc[i]*Pn[id1];

                plotVortx[cnt] = (1-cc[i])*plotVortx[id0] + cc[i]*plotVortx[id1];
                plotVorty[cnt] = (1-cc[i])*plotVorty[id0] + cc[i]*plotVorty[id1];
                plotVortz[cnt] = (1-cc[i])*plotVortz[id0] + cc[i]*plotVortz[id1];
                plotDiv[cnt] = (1-cc[i])*plotDiv[id0] + cc[i]*plotDiv[id1];
                cnt++;
              }
            }
          } else if ((ccnt[6]==4)){
            elementType[NslicedElements++] = 3;
            int m=0;
            for (int i=0;i<6;i++) {
              if (ccnt[i]) {
                int id0 = cid0[i], id1 = cid1[i];
                plotx[cnt] = (1-cc[i])*xn[id0] + cc[i]*xn[id1];
                ploty[cnt] = (1-cc[i])*yn[id0] + cc[i]*yn[id1];
                plotz[cnt] = (1-cc[i])*zn[id0] + cc[i]*zn[id1];
                
                plotU[cnt] = (1-cc[i])*Un[id0] + cc[i]*Un[id1];
                plotV[cnt] = (1-cc[i])*Vn[id0] + cc[i]*Vn[id1];
                plotW[cnt] = (1-cc[i])*Wn[id0] + cc[i]*Wn[id1];
                plotP[cnt] = (1-cc[i])*Pn[id0] + cc[i]*Pn[id1];

                plotVortx[cnt] = (1-cc[i])*plotVortx[id0] + cc[i]*plotVortx[id1];
                plotVorty[cnt] = (1-cc[i])*plotVorty[id0] + cc[i]*plotVorty[id1];
                plotVortz[cnt] = (1-cc[i])*plotVortz[id0] + cc[i]*plotVortz[id1];
                plotDiv[cnt] = (1-cc[i])*plotDiv[id0] + cc[i]*plotDiv[id1];
                cnt++;
                m++;
                if (m==3) break;
              }
            }
            elementType[NslicedElements++] = 3;
            m=0;
            for (int i=0;i<6;i++) {
              if (ccnt[i]) {
                if (m==0) {
                  m++;
                  continue;
                }
                int id0 = cid0[i], id1 = cid1[i];
                plotx[cnt] = (1-cc[i])*xn[id0] + cc[i]*xn[id1];
                ploty[cnt] = (1-cc[i])*yn[id0] + cc[i]*yn[id1];
                plotz[cnt] = (1-cc[i])*zn[id0] + cc[i]*zn[id1];
                
                plotU[cnt] = (1-cc[i])*Un[id0] + cc[i]*Un[id1];
                plotV[cnt] = (1-cc[i])*Vn[id0] + cc[i]*Vn[id1];
                plotW[cnt] = (1-cc[i])*Wn[id0] + cc[i]*Wn[id1];
                plotP[cnt] = (1-cc[i])*Pn[id0] + cc[i]*Pn[id1];

                plotVortx[cnt] = (1-cc[i])*plotVortx[id0] + cc[i]*plotVortx[id1];
                plotVorty[cnt] = (1-cc[i])*plotVorty[id0] + cc[i]*plotVorty[id1];
                plotVortz[cnt] = (1-cc[i])*plotVortz[id0] + cc[i]*plotVortz[id1];
                plotDiv[cnt] = (1-cc[i])*plotDiv[id0] + cc[i]*plotDiv[id1];
                cnt++;
                m++;
              }
            }
          } 
        }
      }
    }

    //write sliced data to file
    FILE *fp;
    fp = fopen(fileName, "w");

    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", cnt, NslicedElements);
    
    // write out nodes
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
    for(iint n=0;n<cnt;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotx[n],ploty[n],plotz[n]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    
    // write out pressure
    fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");
    for(iint n=0;n<cnt;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", plotP[n]);
    }
    fprintf(fp, "       </DataArray>\n");

    // write out vorticity
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"VorticityDivergence\" NumberOfComponents=\"4\" Format=\"ascii\">\n");
    for(iint n=0;n<cnt;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g %g\n", plotVortx[n], plotVorty[n], plotVortz[n], plotDiv[n]);
    }
    fprintf(fp, "       </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
    for(iint n=0;n<cnt;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotU[n], plotV[n], plotW[n]);
    }
    fprintf(fp, "       </DataArray>\n");
    fprintf(fp, "     </PointData>\n");
    
    fprintf(fp, "    <Cells>\n");
    fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
    iint ccnt = 0;
    for(iint n=0;n<NslicedElements;n++){
      fprintf(fp, "       ");
      if (elementType[n]==3) {
        for(int m=0;m<3;++m){
          fprintf(fp, "%d ", ccnt+m);
        }
        ccnt+=3;
      } else {
        for(int m=0;m<4;++m){
          fprintf(fp, "%d ", ccnt+m);
        }
        ccnt+=4;
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");
    
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
    ccnt = 0;
    for(iint n=0;n<NslicedElements;n++){
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", ccnt);
      if (elementType[n]==3) ccnt += 3;
      if (elementType[n]==4) ccnt += 4;
    }
    fprintf(fp, "       </DataArray>\n");
    
    fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
    for(iint n=0;n<NslicedElements;++n){
      if (elementType[n]==3) 
        fprintf(fp, "5\n"); 
      else 
        fprintf(fp, "9\n"); 
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);


    free(elementType);
    free(plotx);
    free(ploty);
    free(plotz);
    free(xn);
    free(yn);
    free(zn);

    free(plotU);
    free(plotV);
    free(plotW);
    free(plotP);
    free(Un);
    free(Vn);
    free(Wn);
    free(Pn);

    free(plotVortx);
    free(plotVorty);
    free(plotVortz);
    free(plotDiv);
    free(Vortx);
    free(Vorty);
    free(Vortz);
    free(Div);
    free(Vortxn);
    free(Vortyn);
    free(Vortzn);
    free(Divn);
#endif

  free(sliceFlag);
  free(sliceSubFlag);
}
