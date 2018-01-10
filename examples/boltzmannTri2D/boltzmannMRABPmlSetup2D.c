#include "boltzmann2D.h"

void boltzmannMRABPmlSetup2D(mesh2D *mesh, char *options){

  //constant pml absorption coefficient
  dfloat xsigma = 80,  ysigma  = 80;
  dfloat cxsigma = 1000, cysigma = 1000;

  //construct element and halo lists
  mesh->MRABpmlNelements = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));
  mesh->MRABpmlElementIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABpmlIds = (iint **) calloc(mesh->MRABNlevels, sizeof(iint*));

  mesh->MRABpmlNhaloElements = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));
  mesh->MRABpmlHaloElementIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABpmlHaloIds = (iint **) calloc(mesh->MRABNlevels, sizeof(iint*));

  //count the pml elements
  mesh->pmlNelements=0;
  for (iint lev =0;lev<mesh->MRABNlevels;lev++){
    for (iint m=0;m<mesh->MRABNelements[lev];m++) {
      iint e = mesh->MRABelementIds[lev][m];
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)) {
        mesh->pmlNelements++;
        mesh->MRABpmlNelements[lev]++;
      }
    }
    for (iint m=0;m<mesh->MRABNhaloElements[lev];m++) {
      iint e = mesh->MRABhaloIds[lev][m];
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300))
        mesh->MRABpmlNhaloElements[lev]++;
    }
  }


  //set up the pml
  if (mesh->pmlNelements) {

    //construct a numbering of the pml elements
    iint *pmlIds = (iint *) calloc(mesh->Nelements,sizeof(iint));
    iint pmlcnt = 0;
    for (iint e=0;e<mesh->Nelements;e++) {
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300))  //pml element
        pmlIds[e] = pmlcnt++;
    }

    //set up lists of pml elements and remove the pml elements from the nonpml MRAB lists
    for (iint lev =0;lev<mesh->MRABNlevels;lev++){
      mesh->MRABpmlElementIds[lev] = (iint *) calloc(mesh->MRABpmlNelements[lev],sizeof(iint));
      mesh->MRABpmlIds[lev] = (iint *) calloc(mesh->MRABpmlNelements[lev],sizeof(iint));
      mesh->MRABpmlHaloElementIds[lev] = (iint *) calloc(mesh->MRABpmlNhaloElements[lev],sizeof(iint));
      mesh->MRABpmlHaloIds[lev] = (iint *) calloc(mesh->MRABpmlNhaloElements[lev],sizeof(iint));

      iint pmlcnt = 0;
      iint nonpmlcnt = 0;
      for (iint m=0;m<mesh->MRABNelements[lev];m++){
        iint e = mesh->MRABelementIds[lev][m];
        int type = mesh->elementInfo[e];

        if ((type==100)||(type==200)||(type==300)) { //pml element
          mesh->MRABpmlElementIds[lev][pmlcnt] = e;
          mesh->MRABpmlIds[lev][pmlcnt] = pmlIds[e];
          pmlcnt++;
        } else { //nonpml element
          mesh->MRABelementIds[lev][nonpmlcnt] = e;
          nonpmlcnt++;
        }
      }

      pmlcnt = 0;
      nonpmlcnt = 0;
      for (iint m=0;m<mesh->MRABNhaloElements[lev];m++){
        iint e = mesh->MRABhaloIds[lev][m];
        int type = mesh->elementInfo[e];

        if ((type==100)||(type==200)||(type==300)) { //pml element
          mesh->MRABpmlHaloElementIds[lev][pmlcnt] = e;
          mesh->MRABpmlHaloIds[lev][pmlcnt] = pmlIds[e];
          pmlcnt++;
        } else { //nonpml element
          mesh->MRABhaloIds[lev][nonpmlcnt] = e;
          nonpmlcnt++;
        }
      }

      //resize nonpml element lists
      mesh->MRABNelements[lev] -= mesh->MRABpmlNelements[lev];
      mesh->MRABNhaloElements[lev] -= mesh->MRABpmlNhaloElements[lev];
      mesh->MRABelementIds[lev] = (iint*) realloc(mesh->MRABelementIds[lev],mesh->MRABNelements[lev]*sizeof(iint));
      mesh->MRABhaloIds[lev]    = (iint*) realloc(mesh->MRABhaloIds[lev],mesh->MRABNhaloElements[lev]*sizeof(iint));
    }

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
    for (iint lev =0;lev<mesh->MRABNlevels;lev++){
      for (iint m=0;m<mesh->MRABpmlNelements[lev];m++) {
        iint e = mesh->MRABpmlElementIds[lev][m];
        iint pmlId = mesh->MRABpmlIds[lev][m];
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
        }
        }
      }
    }

    printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
    mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

    mesh->pmlNfields = 6;

      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsqy);


    mesh->o_pmlSigmaX = mesh->device.malloc(mesh->pmlNelements*mesh->Np*sizeof(dfloat),mesh->pmlSigmaX);
    mesh->o_pmlSigmaY = mesh->device.malloc(mesh->pmlNelements*mesh->Np*sizeof(dfloat),mesh->pmlSigmaY);



    mesh->o_MRABpmlElementIds     = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlIds            = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloElementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloIds        = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
      if (mesh->MRABpmlNelements[lev]) {
        mesh->o_MRABpmlElementIds[lev] = mesh->device.malloc(mesh->MRABpmlNelements[lev]*sizeof(iint),
           mesh->MRABpmlElementIds[lev]);
        mesh->o_MRABpmlIds[lev] = mesh->device.malloc(mesh->MRABpmlNelements[lev]*sizeof(iint),
           mesh->MRABpmlIds[lev]);
      }
      if (mesh->MRABpmlNhaloElements[lev]) {
        mesh->o_MRABpmlHaloElementIds[lev] = mesh->device.malloc(mesh->MRABpmlNhaloElements[lev]*sizeof(iint),
           mesh->MRABpmlHaloElementIds[lev]);
        mesh->o_MRABpmlHaloIds[lev] = mesh->device.malloc(mesh->MRABpmlNhaloElements[lev]*sizeof(iint),
           mesh->MRABpmlHaloIds[lev]);
      }
    }

    free(pmlIds);
  }
}
