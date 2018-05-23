#include "bns.h"

void bnsPmlSetup(bns_t *bns, setupAide &options){

  mesh_t *mesh = bns->mesh;  
  //count the pml elements
  mesh->pmlNelements=0;

  if(options.compareArgs("ABSORBING LAYER", "PML")){
  for (dlong m=0;m<mesh->Nelements;m++) {
    dlong e    = mesh->nonPmlElementIds[m];
    int type   = mesh->elementInfo[e];
    if ((type==100)||(type==200)||(type==300)|| 
        (type==400)||(type==500)||(type==600)||
        (type==700) ) {
      mesh->pmlNelements++;
    }
  }
}

  //set up the pml
  if (mesh->pmlNelements) {

    //construct a numbering of the pml elements
    dlong *pmlIds   = (dlong *) calloc(mesh->Nelements,sizeof(dlong));
    dlong pmlcnt    = 0;
    dlong nonpmlcnt = 0;
    for (dlong e=0;e<mesh->Nelements;e++) {
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)|| 
          (type==400)||(type==500)||(type==600)||
          (type==700) )  //pml element
        pmlIds[e] = pmlcnt++;
    }

    mesh->pmlElementIds = (dlong *) calloc(mesh->pmlNelements,sizeof(dlong));
    mesh->pmlIds        = (dlong *) calloc(mesh->pmlNelements,sizeof(dlong));

    pmlcnt = 0;
    nonpmlcnt = 0;
    for (dlong m=0;m<mesh->Nelements;m++){
      dlong e = mesh->nonPmlElementIds[m];
      int type = mesh->elementInfo[e];

      if ((type==100)||(type==200)||(type==300)|| 
          (type==400)||(type==500)||(type==600)||
          (type==700) ) { //pml element
        mesh->pmlElementIds[pmlcnt] = e;
        mesh->pmlIds[pmlcnt] = pmlIds[e];
        pmlcnt++;
      } else { //nonpml element
        mesh->nonPmlElementIds[nonpmlcnt] = e;
        nonpmlcnt++;
      }
    }


    //resize nonpml element lists
    mesh->nonPmlNelements -= mesh->pmlNelements;
    mesh->nonPmlElementIds = (dlong*) realloc(mesh->nonPmlElementIds,mesh->nonPmlNelements*sizeof(dlong));


   
    //
    printf("Setting PML Coefficient for Cubature Integration\n");
    //set up damping parameter
    bns->pmlSigmaX   = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat));
    bns->pmlSigmaY   = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat)); 
    bns->pmlSigmaZ   = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat)); 

       
        
    //find the bounding box of the whole domain and interior domain
    dfloat xmin = 1e9, xmax =-1e9;
    dfloat ymin = 1e9, ymax =-1e9;
    dfloat zmin = 1e9, zmax =-1e9;
    dfloat pmlxmin = 1e9, pmlxmax =-1e9;
    dfloat pmlymin = 1e9, pmlymax =-1e9;
    dfloat pmlzmin = 1e9, pmlzmax =-1e9;

    for (dlong e=0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Nverts;n++) {
        dfloat x =0, y=0, z = 0; 
        x  = mesh->EX[e*mesh->Nverts+n];
        y  = mesh->EY[e*mesh->Nverts+n];
        
        pmlxmin = (pmlxmin > x) ? x : pmlxmin;
        pmlymin = (pmlymin > y) ? y : pmlymin;
        pmlxmax = (pmlxmax < x) ? x : pmlxmax;
        pmlymax = (pmlymax < y) ? y : pmlymax;

        if(bns->dim==3){
          z = mesh->EZ[e*mesh->Nverts+n];
          pmlzmin = (pmlzmin > z) ? z : pmlzmin;
          pmlzmax = (pmlzmax < z) ? z : pmlzmax;
        }
      }
      //skip pml elements
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)|| 
          (type==400)||(type==500)||(type==600)||
          (type==700) ) continue;
      for (int n=0;n<mesh->Nverts;n++) {
        dfloat x = 0., y = 0., z = 0.; 
        x = mesh->EX[e*mesh->Nverts+n];
        y = mesh->EY[e*mesh->Nverts+n];
        xmin = (xmin > x) ? x : xmin;
        ymin = (ymin > y) ? y : ymin;
        xmax = (xmax < x) ? x : xmax;
        ymax = (ymax < y) ? y : ymax;
        if(bns->dim==3){
          z = mesh->EZ[e*mesh->Nverts+n];
          zmin = (zmin > z) ? z : zmin;
          zmax = (zmax < z) ? z : zmax;          
        }

      }
    }

    dfloat xmaxScale =0., xminScale=0.; 
    dfloat ymaxScale =0., yminScale=0.; 
    dfloat zmaxScale =0., zminScale=0.; 

    xmaxScale = pow(pmlxmax-xmax,bns->pmlOrder);
    xminScale = pow(pmlxmin-xmin,bns->pmlOrder);
    ymaxScale = pow(pmlymax-ymax,bns->pmlOrder);
    yminScale = pow(pmlymin-ymin,bns->pmlOrder);
    if(bns->dim==3){
      zmaxScale = pow(pmlzmax-zmax,bns->pmlOrder);
      zminScale = pow(pmlzmin-zmin,bns->pmlOrder);      
    }

    dfloat *xe  = (dfloat *) calloc(mesh->Nverts, sizeof(dfloat));
    dfloat *ye  = (dfloat *) calloc(mesh->Nverts, sizeof(dfloat));
    dfloat *ze  = (dfloat *) calloc(mesh->Nverts, sizeof(dfloat));

    //set up the damping factor
    for (dlong es=0;es<mesh->pmlNelements;es++){
        dlong e     = mesh->pmlElementIds[es];
        dlong pmlId = mesh->pmlIds[es];
        int type    = mesh->elementInfo[e];
        dlong id    = e*mesh->Nverts;

        for(int n=0; n<mesh->Nverts; n++){
          xe[n] = mesh->EX[id+n];
          ye[n] = mesh->EY[id+n];
          if(bns->dim==3)
            ze[n] = mesh->EZ[id+n];
        }

        for(int n=0;n<mesh->cubNp;++n){ /* for each node */
          dfloat x  = 0, y  = 0, z  = 0; 
          dfloat rn = 0, sn = 0, tn = 0;
          if(bns->elementType==TRIANGLES){
            rn = mesh->cubr[n];
            sn = mesh->cubs[n]; 
            x = -0.5*(rn+sn)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2];
            y = -0.5*(rn+sn)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2];
          }
          else if(bns->elementType==QUADRILATERALS){
            // r is fastest
            const int i = n%mesh->cubNq;
            const int j = n/mesh->cubNq; 
            rn = mesh->cubr[i];
            sn = mesh->cubr[j]; 
            // sn = mesh->cubs[n]; 
            x =  0.25*( (1.0-rn)*(1-sn)*xe[0]+(1.0-rn)*(1+sn)*xe[1]+(1.0+rn)*(1+sn)*xe[2]+(1.0+rn)*(1-sn)*xe[3]);
            y =  0.25*( (1.0-rn)*(1-sn)*ye[0]+(1.0-rn)*(1+sn)*ye[1]+(1.0+rn)*(1+sn)*ye[2]+(1.0+rn)*(1-sn)*ye[3]); 
          }
          else if(bns->elementType==TETRAHEDRA)
          {
            rn = mesh->cubr[n];
            sn = mesh->cubs[n]; 
            tn = mesh->cubt[n];
            x = -0.5*(rn+sn+tn+1)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2] + 0.5*(tn+1)*xe[3];
            y = -0.5*(rn+sn+tn+1)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2] + 0.5*(tn+1)*ye[3];
            z = -0.5*(rn+sn+tn+1)*ze[0] + 0.5*(1+rn)*ze[1] + 0.5*(1+sn)*ze[2] + 0.5*(tn+1)*ze[3];
          } 
                      
          if (type==100) { //X Pml
            if(x>xmax)
              bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmax,bns->pmlOrder)/xmaxScale;
            if(x<xmin)
              bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmin,bns->pmlOrder)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymax,bns->pmlOrder)/ymaxScale;
            if(y<ymin)
              bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymin,bns->pmlOrder)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmax,bns->pmlOrder)/xmaxScale;
            if(x<xmin)
              bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmin,bns->pmlOrder)/xminScale;
            if(y>ymax)
              bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymax,bns->pmlOrder)/ymaxScale;
            if(y<ymin)
              bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymin,bns->pmlOrder)/yminScale;
          }

          if(bns->dim==3){
            if (type==400) { //Z Pml
              if(z>zmax)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmax,bns->pmlOrder)/zmaxScale;
              if(z<zmin)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmin,bns->pmlOrder)/zminScale;
            } else if (type==500) {//XZ Pml
              if(x>xmax)
                bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmax,bns->pmlOrder)/xmaxScale;
              if(x<xmin)
                bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmin,bns->pmlOrder)/xminScale;
              if(z>zmax)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmax,bns->pmlOrder)/zmaxScale;
              if(z<zmin)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmin,bns->pmlOrder)/zminScale;
            } else if (type==600){ //YZ Pml
              if(y>ymax)
                bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymax,bns->pmlOrder)/ymaxScale;
              if(y<ymin)
                bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymin,bns->pmlOrder)/yminScale;
              if(z>zmax)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmax,bns->pmlOrder)/zmaxScale;
              if(z<zmin)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmin,bns->pmlOrder)/zminScale;
            } else if (type==700){ //XYZ Pml
              if(x>xmax)
                bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmax,bns->pmlOrder)/xmaxScale;
              if(x<xmin)
                bns->pmlSigmaX[mesh->cubNp*pmlId + n] = bns->sigmaXmax*pow(x-xmin,bns->pmlOrder)/xminScale;
              if(y>ymax)
                bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymax,bns->pmlOrder)/ymaxScale;
              if(y<ymin)
                bns->pmlSigmaY[mesh->cubNp*pmlId + n] = bns->sigmaYmax*pow(y-ymin,bns->pmlOrder)/yminScale;
              if(z>zmax)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmax,bns->pmlOrder)/zmaxScale;
              if(z<zmin)
                bns->pmlSigmaZ[mesh->cubNp*pmlId + n] = bns->sigmaZmax*pow(z-zmin,bns->pmlOrder)/zminScale;
            }
          }
   
      }

    }

    free(xe); free(ye); free(ze);

    printf("# of PML elements: %d and # of Non-PML elements: %d \n",mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

   

    if(options.compareArgs("TIME INTEGRATOR","LSERK")){
      bns->pmlqx     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlresqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqy     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlresqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqz     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqz  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlresqz  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);
      bns->o_pmlresqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlresqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
      bns->o_pmlresqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlresqy);
      //
      // set up PML on DEVICE    
      bns->o_pmlqz     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqz);
      bns->o_pmlrhsqz  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqz);
      bns->o_pmlresqz  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlresqz);
    }


    if(options.compareArgs("TIME INTEGRATOR","SARK") ){
      bns->pmlqx     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      
      bns->pmlqy     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqz     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqz  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkqx      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkqy      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkqz      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkrhsqx   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkrhsqy   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkrhsqz   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    
      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_saveqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_saveqy    = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);

      bns->o_pmlqz     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqz);
      bns->o_saveqz    = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqz);
      bns->o_pmlrhsqz  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqz);
      
      bns->o_rkqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqx);
      bns->o_rkqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqy);
      bns->o_rkqz     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqz);
      
      bns->o_rkrhsqx  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqx);
      bns->o_rkrhsqy  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqy);
      bns->o_rkrhsqz  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqz);

    }

    if (mesh->pmlNelements){
      bns->o_pmlSigmaX     = mesh->device.malloc(mesh->pmlNelements*mesh->cubNp*sizeof(dfloat),bns->pmlSigmaX);
      bns->o_pmlSigmaY     = mesh->device.malloc(mesh->pmlNelements*mesh->cubNp*sizeof(dfloat),bns->pmlSigmaY);
      bns->o_pmlSigmaZ     = mesh->device.malloc(mesh->pmlNelements*mesh->cubNp*sizeof(dfloat),bns->pmlSigmaZ);

      mesh->o_pmlElementIds = mesh->device.malloc(mesh->pmlNelements*sizeof(dlong), mesh->pmlElementIds);
      mesh->o_pmlIds        = mesh->device.malloc(mesh->pmlNelements*sizeof(dlong), mesh->pmlIds);
    }


      
   free(pmlIds);



  #if 0
       //
 //    printf("Wwitting Pml profile\n" );
 //     int fid = 0; 
 // for (int es=0;es<mesh->pmlNelements;es++){
 //       int e     = mesh->pmlElementIds[es];
 //       int pmlId = mesh->pmlIds[es];   
 //    for(int n=0;n<mesh->Np;n++){
 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 0] = mesh->pmlBetaX[mesh->Np*e+n];
 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 1] = mesh->pmlBetaY[mesh->Np*e+n];
 //    }
 //  }
 
  char fname[BUFSIZ];
  sprintf(fname, "pmlProfile1.vtu");
  boltzmannPlotVTU2D(mesh, fname);

  //    //
 //    printf("Wwitting Pml profile\n" );
 //     int fid = 0; 
 // for (int es=0;es<mesh->pmlNelements;es++){
 //       int e     = mesh->pmlElementIds[es];
 //       int pmlId = mesh->pmlIds[es];   
 //    for(int n=0;n<mesh->cubNp;n++){
 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 0] = bns->pmlSigmaX[mesh->Np*e+n];
 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 1] = bns->pmlSigmaY[mesh->Np*e+n];
 //    }
 //  }


 //  char fname[BUFSIZ];
 //  sprintf(fname, "pmlProfile1.vtu");
 //  boltzmannPlotVTU2D(mesh, fname);



 //   dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
 //   mesh->o_cubProjectT.copyTo(cubProjectT);
 //    // free(pmlIds);
 //     int fid = 0; 
 // for (int es=0;es<mesh->pmlNelements;es++){
 //       int e     = mesh->pmlElementIds[es];
 //       int pmlId = mesh->pmlIds[es];   
 //    for(int n=0;n<mesh->Np;n++)
 //    {
 //      dfloat q1 = 0; dfloat q2 = 0; 
 //      for(int m=0; m<mesh->cubNp;++m){
 //         dfloat prj = cubProjectT[m*mesh->Np+n];
 //         q1 += prj*bns->pmlSigmaX[mesh->cubNp*pmlId + m];
 //         q2 += prj*bns->pmlSigmaY[mesh->cubNp*pmlId + m];
 //      }

 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 0] = q1;
 //      mesh->q[bns->Nfields*(n + e*mesh->Np) + 1] = q2;
 //    }
 //  }
  

 //  // char fname[BUFSIZ];
 //  // sprintf(fname, "pmlProfile1.vtu");
 //  // boltzmannPlotVTU2D(mesh, fname);

 //   //
 //  printf("NonPmlElements: %d  PmlElements: %d \n", mesh->nonPmlNelements, mesh->pmlNelements);

  
 //  for (int es=0;es<mesh->nonPmlNelements;es++){
 //       int e      = mesh->nonPmlElementIds[es];
 //       int id = e*mesh->Nverts;  
 //        dfloat x1 = mesh->EX[id+0]; /* x-coordinates of vertices */
 //        dfloat x2 = mesh->EX[id+1];
 //        dfloat x3 = mesh->EX[id+2];

 //        dfloat y1 = mesh->EY[id+0]; /* y-coordinates of vertices */
 //        dfloat y2 = mesh->EY[id+1];
 //        dfloat y3 = mesh->EY[id+2];
 //       printf("Element: %d  x1: %.5f x2: %.5f x3: %.5f y1: %.5f y2: %.5f y3: %.5f \n", e, x1,x2,x3,y1,y2,y3);
 //  }   


 //  for (int es=0;es<mesh->pmlNelements;es++){
 //       int e     = mesh->pmlElementIds[es];
 //       int pmlId = mesh->pmlIds[es];
 //       int id = e*mesh->Nverts;  
 //       //
 //       dfloat x1 = mesh->EX[id+0]; /* x-coordinates of vertices */
 //        dfloat x2 = mesh->EX[id+1];
 //        dfloat x3 = mesh->EX[id+2];

 //        dfloat y1 = mesh->EY[id+0]; /* y-coordinates of vertices */
 //        dfloat y2 = mesh->EY[id+1];
 //        dfloat y3 = mesh->EY[id+2];
 //       printf("Element: %d   PmlId = %d x1: %.5f x2: %.5f x3: %.5f y1: %.5f y2: %.5f y3: %.5f\n", e, pmlId, x1,x2,x3,y1,y2,y3);
 //  }   



  #endif
  
  }

}
