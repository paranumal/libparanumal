#include "acoustics2D.h"

void acousticsPmlSetup2D(mesh2D *mesh){

  //constant pml absorption coefficient
  dfloat xsigma = 80, ysigma = 80;

  //construct element and halo lists
  mesh->MRABpmlNelements     = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));
  mesh->MRABpmlNhaloElements = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));
  
  mesh->MRABpmlElementIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABpmlIds        = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));

  mesh->MRABpmlHaloElementIds = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABpmlHaloIds        = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));

  mesh->MRABpmlNelP      = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  mesh->MRABpmlNhaloEleP = (iint **) calloc(mesh->MRABNlevels,sizeof(iint*));
  
  mesh->MRABpmlElIdsP = (iint ***) calloc(mesh->MRABNlevels,sizeof(iint**));
  mesh->MRABpmlIdsP   = (iint ***) calloc(mesh->MRABNlevels,sizeof(iint**));  
  
  mesh->MRABpmlHaloEleIdsP = (iint ***) calloc(mesh->MRABNlevels,sizeof(iint**));
  mesh->MRABpmlHaloIdsP    = (iint ***) calloc(mesh->MRABNlevels,sizeof(iint**));
  
  for (iint lev=0;lev<mesh->MRABNlevels;lev++) {
    mesh->MRABpmlNelP[lev]      = (iint *) calloc(mesh->NMax+1,sizeof(iint));
    mesh->MRABpmlNhaloEleP[lev] = (iint *) calloc(mesh->NMax+1,sizeof(iint));
    mesh->MRABpmlElIdsP = (iint **) calloc(mesh->NMax+1,sizeof(iint *));
    mesh->MRABpmlIdsP   = (iint **) calloc(mesh->NMax+1,sizeof(iint *));  
    mesh->MRABpmlHaloEleIdsP = (iint **) calloc(mesh->NMax+1,sizeof(iint *));
    mesh->MRABpmlHaloIdsP    = (iint **) calloc(mesh->NMax+1,sizeof(iint *));
  }  

  //count the pml elements
  mesh->pmlNelements=0;
  for (iint lev =0;lev<mesh->MRABNlevels;lev++){
    for (iint m=0;m<mesh->MRABNelements[lev];m++) {
      iint e = mesh->MRABelementIds[lev][m];
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)) {
        mesh->pmlNelements++;
        mesh->MRABpmlNelements[lev]++;
        mesh->MRABpmlNelP[lev][mesh->N[e]]++;
      }
    }
    for (iint m=0;m<mesh->MRABNhaloElements[lev];m++) {
      iint e = mesh->MRABhaloIds[lev][m];
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300)) {
        mesh->MRABpmlNhaloElements[lev]++;
        mesh->MRABpmlNhaloEleP[lev][mesh->N[e]]++;
      }
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
      mesh->MRABpmlIds[lev]        = (iint *) calloc(mesh->MRABpmlNelements[lev],sizeof(iint));
      mesh->MRABpmlHaloElementIds[lev] = (iint *) calloc(mesh->MRABpmlNhaloElements[lev],sizeof(iint));
      mesh->MRABpmlHaloIds[lev]        = (iint *) calloc(mesh->MRABpmlNhaloElements[lev],sizeof(iint));

      for (iint p=0;p<=mesh->NMax;p++) {
        mesh->MRABpmlElIdsP[lev][p] = (iint *) calloc(mesh->MRABpmlNelP[lev][p],sizeof(iint));
        mesh->MRABpmlIdsP[lev][p]   = (iint *) calloc(mesh->MRABpmlNelP[lev][p],sizeof(iint));
        mesh->MRABpmlHaloEleIdsP[lev][p] = (iint *) calloc(mesh->MRABpmlNhaloEleP[lev][p],sizeof(iint)); 
        mesh->MRABpmlHaloIdsP[lev][p]    = (iint *) calloc(mesh->MRABpmlNhaloEleP[lev][p],sizeof(iint)); 
      }

      //initialize counters
      iint pmlcnt = 0;
      iint nonpmlcnt = 0;
      iint pmlPcnt[mesh->NMax+1];
      iint nonpmlPcnt[mesh->NMax+1];
      for (iint p=0;p<=mesh->NMax;p++) {
        pmlPcnt[p] = 0;
        nonpmlPcnt[p] = 0;
      }
      for (iint m=0;m<mesh->MRABNelements[lev];m++){
        iint e = mesh->MRABelementIds[lev][m];
        int N = mesh->N[e];
        int type = mesh->elementInfo[e];

        if ((type==100)||(type==200)||(type==300)) { //pml element
          mesh->MRABpmlElementIds[lev][pmlcnt] = e;
          mesh->MRABpmlIds[lev][pmlcnt] = pmlIds[e];
          pmlcnt++;
          mesh->MRABpmlEleIdsP[lev][N][pmlPcnt[N]] = e;
          mesh->MRABpmlIdsP[lev][N][pmlPcnt[N]] = pmlIds[e];
          pmlPcnt[N]++;
        } else { //nonpml element
          mesh->MRABelementIds[lev][nonpmlcnt] = e;
          nonpmlcnt++;
          mesh->MRABelIdsP[lev][N][nonpmlcnt[N]] = e;
          nonpmlPcnt[N]++;
        }
      }

      pmlcnt = 0;
      nonpmlcnt = 0;
      for (iint p=0;p<=mesh->NMax;p++) {
        pmlPcnt[p] = 0;
        nonpmlPcnt[p] = 0;
      }
      for (iint m=0;m<mesh->MRABNhaloElements[lev];m++){
        iint e = mesh->MRABhaloIds[lev][m];
        int N = mesh->N[e];
        int type = mesh->elementInfo[e];

        if ((type==100)||(type==200)||(type==300)) { //pml element
          mesh->MRABpmlHaloElementIds[lev][pmlcnt] = e;
          mesh->MRABpmlHaloIds[lev][pmlcnt] = pmlIds[e];
          pmlcnt++;
          mesh->MRABpmlHaloEleIdsP[lev][N][pmlPcnt[N]] = e;
          mesh->MRABpmlHaloIdsP[lev][N][pmlPcnt[N]] = pmlIds[e];
          pmlPcnt[N]++;
        } else { //nonpml element
          mesh->MRABhaloIds[lev][nonpmlcnt] = e;
          nonpmlcnt++;
          mesh->MRABhaloIdsP[lev][N][nonpmlPcnt[N]] = e;
          nonpmlPcnt[N]++;
        }
      }

      //resize nonpml element lists
      mesh->MRABNelements[lev]     -= mesh->MRABpmlNelements[lev];
      mesh->MRABNhaloElements[lev] -= mesh->MRABpmlNhaloElements[lev];
      mesh->MRABelementIds[lev] = (iint*) realloc(mesh->MRABelementIds[lev],mesh->MRABNelements[lev]*sizeof(iint));
      mesh->MRABhaloIds[lev]    = (iint*) realloc(mesh->MRABhaloIds[lev],mesh->MRABNhaloElements[lev]*sizeof(iint));
    
      for (iint p=0;p<=mesh->NMax;p++) {
        mesh->MRABNelP[lev][p]      -= mesh->MRABpmlNelP[lev][p];
        mesh->MRABNhaloEleP[lev][p] -= mesh->MRABpmlNhaloEleP[lev][p];
        mesh->MRABelIdsP[lev][p]   = (iint*) realloc(mesh->MRABelIdsP[lev][p],mesh->MRABNelP[lev][p]*sizeof(iint));
        mesh->MRABhaloIdsP[lev][p] = (iint*) realloc(mesh->MRABhaloIdsP[lev][p],mesh->MRABNhaloEleP[lev][p]*sizeof(iint));
      }
    }

    //set up damping parameter
    mesh->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNpMax,sizeof(dfloat));
    mesh->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNpMax,sizeof(dfloat));

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
        int N = mesh->N[e];
        iint pmlId = mesh->MRABpmlIds[lev][m];
        int type = mesh->elementInfo[e];

        iint id = e*mesh->Nverts;

        dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = mesh->EX[id+1];
        dfloat xe3 = mesh->EX[id+2];

        dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = mesh->EY[id+1];
        dfloat ye3 = mesh->EY[id+2];

        for(iint n=0;n<mesh->cubNp[N];++n){ /* for each node */
          // cubature node coordinates
          dfloat rn = mesh->cubr[N][n];
          dfloat sn = mesh->cubs[N][n];

          /* physical coordinate of interpolation node */
          dfloat x = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
          dfloat y = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;

          if (type==100) { //X Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->cubNpMax*pmlId + n] = xsigma*pow(x-xmax,2)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->cubNpMax*pmlId + n] = xsigma*pow(x-xmin,2)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlSigmaY[mesh->cubNpMax*pmlId + n] = ysigma*pow(y-ymax,2)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->cubNpMax*pmlId + n] = ysigma*pow(y-ymin,2)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              mesh->pmlSigmaX[mesh->cubNpMax*pmlId + n] = xsigma*pow(x-xmax,2)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[mesh->cubNpMax*pmlId + n] = xsigma*pow(x-xmin,2)/xminScale;
            if(y>ymax)
              mesh->pmlSigmaY[mesh->cubNpMax*pmlId + n] = ysigma*pow(y-ymax,2)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[mesh->cubNpMax*pmlId + n] = ysigma*pow(y-ymin,2)/yminScale;
          }
        }
      }
    }

    printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
    mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

    mesh->pmlNfields = 2;
    mesh->pmlq    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
    mesh->pmlrhsq = (dfloat*) calloc(3*mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

    // set up PML on DEVICE
    mesh->o_pmlq      = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlq);
    mesh->o_pmlrhsq   = mesh->device.malloc(3*mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsq);
    mesh->o_pmlSigmaX = mesh->device.malloc(mesh->pmlNelements*mesh->cubNp*sizeof(dfloat),mesh->pmlSigmaX);
    mesh->o_pmlSigmaY = mesh->device.malloc(mesh->pmlNelements*mesh->cubNp*sizeof(dfloat),mesh->pmlSigmaY);

    mesh->o_MRABpmlElementIds     = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlIds            = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloElementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloIds        = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));

    mesh->o_MRABpmlElIdsP      = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlIdsP        = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloEleIdsP = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    mesh->o_MRABpmlHaloIdsP    = (occa::memory **) malloc(mesh->MRABNlevels*sizeof(occa::memory));

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

      mesh->o_MRABpmlElIdsP[lev]      = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
      mesh->o_MRABpmlIdsP[lev]        = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
      mesh->o_MRABpmlHaloEleIdsP[lev] = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));
      mesh->o_MRABpmlHaloIdsP[lev]    = (occa::memory *) malloc((mesh->NMax+1)*sizeof(occa::memory));        
      for (iint p=0;p<=mesh->NMax;p++) {
        if (mesh->MRABpmlNelP[lev][p]) {
          mesh->o_MRABpmlElIdsP[lev][p] = mesh->device.malloc(mesh->MRABpmlNelP[lev][p]*sizeof(iint),
             mesh->MRABpmlElIdsP[lev][p]);
          mesh->o_MRABpmlIdsP[lev][p] = mesh->device.malloc(mesh->MRABpmlNelP[lev][p]*sizeof(iint),
             mesh->MRABpmlIdsP[lev][p]);
        }
        if (mesh->MRABpmlNhaloEleP[lev][p]) {
          mesh->o_MRABpmlHaloEleIdsP[lev][p] = mesh->device.malloc(mesh->MRABpmlNhaloEleP[lev][p]*sizeof(iint),
             mesh->MRABpmlHaloEleIdsP[lev][p]);
          mesh->o_MRABpmlHaloIdsP[lev][p] = mesh->device.malloc(mesh->MRABpmlNhaloEleP[lev][p]*sizeof(iint),
             mesh->MRABpmlHaloIdsP[lev][p]);
        }
      }
    }

    free(pmlIds);
  }
}
