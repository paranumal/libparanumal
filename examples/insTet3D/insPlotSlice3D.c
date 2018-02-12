#include "ins3D.h"

// interpolate data to plot nodes and save to file (one per process
void insPlotSlice3D(ins_t *ins, char *fileName, const int Nslices, const char** dim, const dfloat* c){

  mesh3D *mesh = ins->mesh;
  
  //find number of sliced elements
  dfloat *coord;
  iint NslicedElements = 0;
  int *sliceFlag = (int *) calloc(mesh->Nelements,sizeof(int));
  
  for (int i=0;i<Nslices;i++) {
    if (strstr(dim[i],"x")) coord = mesh->EX;
    if (strstr(dim[i],"y")) coord = mesh->EY;
    if (strstr(dim[i],"z")) coord = mesh->EZ;

    for (iint e=0;e<mesh->Nelements;e++) {
      iint id = e*mesh->Nverts;
      dfloat coord0 = coord[id+0];
      dfloat coord1 = coord[id+1];
      dfloat coord2 = coord[id+2];
      dfloat coord3 = coord[id+3];

      dfloat c0 = (c[i]-coord0)/(coord1-coord0);
      dfloat c1 = (c[i]-coord1)/(coord2-coord1);
      dfloat c2 = (c[i]-coord2)/(coord0-coord2);
      dfloat c3 = (c[i]-coord0)/(coord3-coord0);
      dfloat c4 = (c[i]-coord1)/(coord3-coord1);
      dfloat c5 = (c[i]-coord2)/(coord3-coord2);

      int cnt =0;
      cnt += ((c0>=0.0)&&(c0<=1.0));
      cnt += ((c1>=0.0)&&(c1<=1.0));
      cnt += ((c2>=0.0)&&(c2<=1.0));
      cnt += ((c3>=0.0)&&(c3<=1.0));
      cnt += ((c4>=0.0)&&(c4<=1.0));
      cnt += ((c5>=0.0)&&(c5<=1.0));

      if (cnt>0) {
        if (sliceFlag[e]==0) NslicedElements++;
        sliceFlag[e] = 1;
      }
    }
  }

  if (NslicedElements) {
    dfloat *Vx   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
    dfloat *Vy   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
    dfloat *Vz   = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
    dfloat *divU = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

    //write sliced data to file
    FILE *fp;
    fp = fopen(fileName, "w");

    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", NslicedElements*mesh->plotNp, NslicedElements*mesh->plotNelements);
    
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

      for(iint n=0;n<mesh->Np;++n){
        dfloat dUdr = 0, dUds = 0, dUdt = 0 ;

        dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
        dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
        dfloat rz = mesh->vgeo[e*mesh->Nvgeo+RZID];    
        
        dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
        dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];
        dfloat sz = mesh->vgeo[e*mesh->Nvgeo+SZID];    
       
        dfloat tx = mesh->vgeo[e*mesh->Nvgeo+TXID];
        dfloat ty = mesh->vgeo[e*mesh->Nvgeo+TYID];
        dfloat tz = mesh->vgeo[e*mesh->Nvgeo+TZID];    

        for(iint m=0;m<mesh->Np;++m){
          iint id = m+e*mesh->Np;
          id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

          dfloat Un = ins->U[id];
          dfloat Vn = ins->V[id];
          dfloat Wn = ins->W[id];

          dUdr += mesh->Dr[n*mesh->Np+m]*(rx*Un+ry*Vn+rz*Wn);
          dUds += mesh->Ds[n*mesh->Np+m]*(sx*Un+sy*Vn+sz*Wn);
          dUdt += mesh->Dt[n*mesh->Np+m]*(tx*Un+ty*Vn+tz*Wn);
        }
        
        // Compute divergence
        divU[n] = dUdr + dUds + dUdt; 
      }

      for(iint n=0;n<mesh->plotNp;++n){
        dfloat plotDiv = 0;
        for(iint m=0;m<mesh->Np;++m){
          plotDiv += mesh->plotInterp[n*mesh->Np+m]*divU[m];
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
        Vx[n] = dWdy-dVdz;
        Vy[n] = dUdz-dWdx;
        Vz[n] = dVdx-dUdy;
      }
      
      for(iint n=0;n<mesh->plotNp;++n){
        dfloat plotVxn = 0, plotVyn = 0, plotVzn = 0 ;
        dfloat plotDivUn = 0;
        for(iint m=0;m<mesh->Np;++m){
          plotVxn   += mesh->plotInterp[n*mesh->Np+m]*Vx[m];
          plotVyn   += mesh->plotInterp[n*mesh->Np+m]*Vy[m];
          plotVzn   += mesh->plotInterp[n*mesh->Np+m]*Vz[m];
          plotDivUn += mesh->plotInterp[n*mesh->Np+m]*divU[m]; 
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
        fprintf(fp, "10\n"); // TET code ?
      }
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);

    free(Vx);
    free(Vy);
    free(Vz);
    free(divU);

#if 0
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

    if (strstr(dim,"x")) coord = xn;
    if (strstr(dim,"y")) coord = yn;
    if (strstr(dim,"z")) coord = zn;

    iint cnt = 0;
    NslicedElements =0;
    for (iint e=0;e<mesh->Nelements;e++) {
      //if (sliceFlag[e]==0) continue;

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
        cid0[3] = id0; cid1[3] = id3;
        cid0[4] = id1; cid1[4] = id3;
        cid0[5] = id2; cid1[5] = id3;

        dfloat cc[6];
        int ccnt =0;
        for (int i=0;i<6;i++) {
          int id0 = cid0[i], id1 = cid1[i];
          cc[i] = (c-coord[id0])/(coord[id1]-coord[id0]);
          ccnt += ((cc[i]>=0.0)&&(cc[i]<=1.0));
        }
        if (ccnt>4) printf("huh, again?\n");
        if ((ccnt==3)||(ccnt==4)) { //if sliced, add the data to the array
          elementType[NslicedElements++] = ccnt;
          for (int i=0;i<6;i++) {
            if ((cc[i]>=0.0)&&(cc[i]<=1.0)) {
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
  }
#endif
  }
  free(sliceFlag);
}
