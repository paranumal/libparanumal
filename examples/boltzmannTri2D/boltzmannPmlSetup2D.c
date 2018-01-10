#include "boltzmann2D.h"

void boltzmannPmlSetup2D(mesh2D *mesh, char *options){

  //constant pml absorption coefficient
  dfloat xsigma = 80.,  ysigma  = 80.;
  dfloat cxsigma =80, cysigma = 80;

  //count the pml elements
  mesh->pmlNelements=0;


    for (iint m=0;m<mesh->Nelements;m++) {
      iint e = mesh->nonPmlElementIds[m];
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)) {
        mesh->pmlNelements++;
      }
    }


printf("Number of pmlElements: %d \n", mesh->pmlNelements);

  //set up the pml
  if (mesh->pmlNelements) {

    //construct a numbering of the pml elements
    iint *pmlIds = (iint *) calloc(mesh->Nelements,sizeof(iint));
    iint pmlcnt  = 0;
    iint nonpmlcnt = 0;
    for (iint e=0;e<mesh->Nelements;e++) {
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300))  //pml element
        pmlIds[e] = pmlcnt++;
    }

    //set up lists of pml elements and remove the pml elements from the nonpml element list

      mesh->pmlElementIds = (iint *) calloc(mesh->pmlNelements,sizeof(iint));
      mesh->pmlIds        = (iint *) calloc(mesh->pmlNelements,sizeof(iint));
     

      pmlcnt = 0;
      nonpmlcnt = 0;
      for (iint m=0;m<mesh->Nelements;m++){
        iint e = mesh->nonPmlElementIds[m];
        int type = mesh->elementInfo[e];

        if ((type==100)||(type==200)||(type==300)) { //pml element
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
      mesh->nonPmlElementIds = (iint*) realloc(mesh->nonPmlElementIds,mesh->nonPmlNelements*sizeof(iint));
     



    //set up damping parameter
    mesh->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));
    mesh->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));

    //find the bounding box of the whole domain and interior domain
    dfloat xmin = 1e9, xmax =-1e9;
    dfloat ymin = 1e9, ymax =-1e9;
    dfloat pmlxmin = 1e9, pmlxmax =-1e9;
    dfloat pmlymin = 1e9, pmlymax =-1e9;
    for (iint e=0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Nverts;n++) {
        dfloat x = mesh->EX[e*mesh->Nverts+n];
        dfloat y = mesh->EY[e*mesh->Nverts+n];

        pmlxmin = (pmlxmin > x) ? x : pmlxmin;
        pmlymin = (pmlymin > y) ? y : pmlymin;
        pmlxmax = (pmlxmax < x) ? x : pmlxmax;
        pmlymax = (pmlymax < y) ? y : pmlymax;
      }

      //skip pml elements
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)) continue;

      for (int n=0;n<mesh->Nverts;n++) {
        dfloat x = mesh->EX[e*mesh->Nverts+n];
        dfloat y = mesh->EY[e*mesh->Nverts+n];

        xmin = (xmin > x) ? x : xmin;
        ymin = (ymin > y) ? y : ymin;
        xmax = (xmax < x) ? x : xmax;
        ymax = (ymax < y) ? y : ymax;
      }
    }






    dfloat xmaxScale = pow(pmlxmax-xmax,2);
    dfloat xminScale = pow(pmlxmin-xmin,2);
    dfloat ymaxScale = pow(pmlymax-ymax,2);
    dfloat yminScale = pow(pmlymin-ymin,2);

    //set up the damping factor
    for (iint es=0;es<mesh->pmlNelements;es++){
        iint e     = mesh->pmlElementIds[es];
        iint pmlId = mesh->pmlIds[es];
        int type = mesh->elementInfo[e];
        iint id = e*mesh->Nverts;

        dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = mesh->EX[id+1];
        dfloat xe3 = mesh->EX[id+2];

        dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = mesh->EY[id+1];
        dfloat ye3 = mesh->EY[id+2];

        for(iint n=0;n<mesh->Np;++n){ /* for each node */
          dfloat x = mesh->x[n + e*mesh->Np];
          dfloat y = mesh->y[n + e*mesh->Np];

         if(strstr(options,"QUADRATIC")){

          if (type==100) { //X Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = xsigma*pow(x-xmax,2)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = xsigma*pow(x-xmin,2)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = ysigma*pow(y-ymax,2)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = ysigma*pow(y-ymin,2)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = xsigma*pow(x-xmax,2)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = xsigma*pow(x-xmin,2)/xminScale;
            if(y>ymax)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = ysigma*pow(y-ymax,2)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = ysigma*pow(y-ymin,2)/yminScale;
          }

        }

        if(strstr(options,"CONSTANT")){

           if (type==100) { //X Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = cxsigma;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = cxsigma;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = cysigma;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = cysigma;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = cxsigma;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->Np*pmlId + n] = cxsigma;
            if(y>ymax)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = cysigma;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->Np*pmlId + n] = cysigma;
          }

         //printf("%.6f %.6f\n", mesh->pmlSigmaX[mesh->Np*pmlId + n],mesh->pmlSigmaY[mesh->Np*pmlId + n]);

        }
        }
      }

    printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
    mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

    mesh->pmlNfields = 6;

    if(strstr(options,"SAAB")  || strstr(options,"SRAB") ){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqy);
    }



    if(strstr(options,"LSERK")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlresqx = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlresqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqx);
      mesh->o_pmlresqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlresqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqy);
      mesh->o_pmlresqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlresqy);
    }


    if(strstr(options,"SARK")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qSx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qSy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqy);
    }


    if(strstr(options,"LSIMEX")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      // set up PML on DEVICE 
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);   
      mesh->o_qSx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qYx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qZx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qSy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qYy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qZy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
    }


    


    mesh->o_pmlSigmaX = mesh->device.malloc(mesh->pmlNelements*mesh->Np*sizeof(dfloat),mesh->pmlSigmaX);
    mesh->o_pmlSigmaY = mesh->device.malloc(mesh->pmlNelements*mesh->Np*sizeof(dfloat),mesh->pmlSigmaY);




      if (mesh->pmlNelements) {
        mesh->o_pmlElementIds = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlElementIds);
        mesh->o_pmlIds        = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlIds);
      }
      
   // free(pmlIds);
  }
}
