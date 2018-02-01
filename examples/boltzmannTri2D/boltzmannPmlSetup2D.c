#include "boltzmann2D.h"

void boltzmannPmlSetup2D(mesh2D *mesh, char *options){

  //constant pml absorption coefficient
  //dfloat xsigma  = 100. , ysigma  = 100.;

  // dfloat Sigma[4]; Sigma[0]=20; Sigma[1] = 50; Sigma[2] = 200; Sigma[3] = 400; 

  // dfloat xsigma  = Sigma[mesh->Ntscale]; dfloat  ysigma  = Sigma[mesh->Ntscale];


  dfloat xsigma  = 200., ysigma  = 200.;
  
  dfloat xsmin = -2.0, xsmax = 2.0; // map x to this range to control erf profile 
  dfloat ysmin = -2.0, ysmax = 2.0; // map y to this range to control erf profile 

  dfloat q1    = 0.2, q2 = 0.8;    // Ramped profile coefficients 0<q1, 1>q2 , q1<exp<q2 or polynomial

// fifth order polynomialcoefficients for q1=0.2 and q2=0.8;
  dfloat c5 =  6.0 ;
  dfloat c4 = -15.0;
  dfloat c3 =  10.0;
  dfloat c2 =  0.0 ;
  dfloat c1 =  0.0 ;
  dfloat c0 =  0.0 ;

  // third order
  // dfloat c3 = -2.0;
  // dfloat c2 =  3.0;
  // dfloat c1 =  0.0;
  // dfloat c0 =  0.0;

  printf(" Sigma Scale: %d\n",mesh->Ntscale);
  //count the pml elements
  mesh->pmlNelements=0;


  for (iint m=0;m<mesh->Nelements;m++) {
    iint e = mesh->nonPmlElementIds[m];
    int type = mesh->elementInfo[e];
    if ((type==100)||(type==200)||(type==300)) {
      mesh->pmlNelements++;
    }
  }

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


    iint Nnodes = 0;
    if(strstr(options,"CUBATURE")){ // !!!!!!!!!!!!!!!
      //
      printf("Setting PML Coefficient for Cubature Integration\n");
      //set up damping parameter
      mesh->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat));
      mesh->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat));  
      Nnodes = mesh->cubNp;    
    }
     else{
      printf("Setting PML Coefficients for Nodal Collocation Integration\n");
      //set up damping parameter
      mesh->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));
      mesh->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));
      Nnodes = mesh->Np; 
    }
   
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



    iint order = 0;
    //
    if(strstr(options,"CONSTANT"))
      order = 0; 
    if(strstr(options,"LINEAR"))
      order = 1; 
    if(strstr(options,"QUADRATIC"))
      order = 2;
    if(strstr(options,"FORTHORDER"))
      order = 4;
    if(strstr(options, "ERFFUNCTION"))
      order =1;

    dfloat xmaxScale = pow(pmlxmax-xmax,order);
    dfloat xminScale = pow(pmlxmin-xmin,order);
    dfloat ymaxScale = pow(pmlymax-ymax,order);
    dfloat yminScale = pow(pmlymin-ymin,order);

    printf("xmaxScale = %.10e \n", xmaxScale);

    //set up the damping factor
    for (iint es=0;es<mesh->pmlNelements;es++){
        iint e     = mesh->pmlElementIds[es];
        iint pmlId = mesh->pmlIds[es];
        iint type  = mesh->elementInfo[e];
        iint id    = e*mesh->Nverts;

        dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = mesh->EX[id+1];
        dfloat xe3 = mesh->EX[id+2];

        dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = mesh->EY[id+1];
        dfloat ye3 = mesh->EY[id+2];

        for(iint n=0;n<Nnodes;++n){ /* for each node */
          dfloat x = 0, y = 0; 

          if(Nnodes==mesh->cubNp){
            dfloat rn = mesh->cubr[n];
            dfloat sn = mesh->cubs[n]; 
            x = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
            y = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
          }
          else{
            x = mesh->x[n + e*mesh->Np];
            y = mesh->y[n + e*mesh->Np];
          }
          
          if(!strstr(options,"ERFUNCTION") && !strstr(options,"RAMPFUNCTION") && !strstr(options,"SMOOTHPOLYNOMIAL")){
          if (type==100) { //X Pml
            if(x>xmax)
              mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmin,order)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymin,order)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmin,order)/xminScale;
            if(y>ymax)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymin,order)/yminScale;
          }
        }

        else if(strstr(options,"ERFUNCTION")){
          if (type==100) { //X Pml
            if(x>xmax)
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmax)/xmaxScale * (xsmax-xsmin)));
            if(x<xmin)
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmin)/xminScale * (xsmax-xsmin)));
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymax)/ymaxScale * (ysmax-ysmin)));
            if(y<ymin)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymin)/yminScale * (ysmax-ysmin)));
          } else if (type==300) { //XY Pml
            if(x>xmax)
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmax)/xmaxScale * (xsmax-xsmin)));
            if(x<xmin)
              mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmin)/xminScale * (xsmax-xsmin)));
            if(y>ymax)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymax)/ymaxScale * (ysmax-ysmin)));
            if(y<ymin)
              mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymin)/yminScale * (ysmax-ysmin)));
          }
        }

        else if(strstr(options,"RAMPFUNCTION")){

          dfloat qx=0,      qy = 0;
          dfloat taux = 0,  tauy = 0;
          dfloat rampx = 0, rampy = 0;
          
          // stretch x to [0 1] interval
          if (type==100) { //X Pml
            if(x>xmax)
              qx   = (x-xmax)/(pmlxmax-xmax);
            if(x<xmin)
              qx    = (x-xmin)/(pmlxmin-xmin);
          } else if (type==200) { //Y Pml
            if(y>ymax)
              qy    = (y-ymax)/(pmlymax-ymax);
            if(y<ymin)
              qy    = (y-ymin)/(pmlymin-ymin);
          } else if (type==300) { //XY Pml
            if(x>xmax)
             qx   = (x-xmax)/(pmlxmax-xmax);
            if(x<xmin)
             qx    = (x-xmin)/(pmlxmin-xmin);
            if(y>ymax)
             qy    = (y-ymax)/(pmlymax-ymax);
            if(y<ymin)
             qy    = (y-ymin)/(pmlymin-ymin);
          }
          
          // Renormalize for ramp profile i.e. tau = [0,1]
          taux  =  (qx-q1)/(q2-q1);
          tauy  =  (qy-q1)/(q2-q1);

          rampx = 1.0 -exp(2.*pow(exp(1.0),(-1. /taux))/(taux-1.0));
          rampy = 1.0 -exp(2.*pow(exp(1.0),(-1. /tauy))/(tauy-1.0));
          rampx = qx<=q1 ? 0. : rampx;  rampx = qx>=q2 ? 1. : rampx;
          rampy = qy<=q1 ? 0. : rampy;  rampy = qy>=q2 ? 1. : rampy;
           
          if (type==100) { //X Pml
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*rampx;
          } else if (type==200) { //Y Pml
             mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*rampy;
          } else if (type==300) { //XY Pml
            if(x>xmax || x<xmin)
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*rampx;
            if(y>ymax || y<ymin)
             mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*rampy;
          }
        }

        else if(strstr(options,"SMOOTHPOLYNOMIAL")){

          dfloat qx=0,      qy = 0;
          dfloat taux = 0,  tauy = 0;
          dfloat polyx = 0, polyy = 0;

          if (type==100) { //X Pml
            if(x>xmax)
              qx   = (x-xmax)/(pmlxmax-xmax);
            if(x<xmin)
              qx    = (x-xmin)/(pmlxmin-xmin);
          } else if (type==200) { //Y Pml
            if(y>ymax)
              qy    = (y-ymax)/(pmlymax-ymax);
            if(y<ymin)
              qy    = (y-ymin)/(pmlymin-ymin);
          } else if (type==300) { //XY Pml
            if(x>xmax)
             qx   = (x-xmax)/(pmlxmax-xmax);
            if(x<xmin)
             qx    = (x-xmin)/(pmlxmin-xmin);
            if(y>ymax)
             qy    = (y-ymax)/(pmlymax-ymax);
            if(y<ymin)
             qy    = (y-ymin)/(pmlymin-ymin);
          }
          // Third Order
          taux  =  (qx-q1)/(q2-q1);
          tauy  =  (qy-q1)/(q2-q1);
          // polyx = c3*pow(taux,3) + c2*pow(taux,2) +c1*taux +c0;
          // polyy = c3*pow(tauy,3) + c2*pow(tauy,2) +c1*tauy +c0;
          // // Fifth Order
          polyx = c5*pow(taux,5) + c4*pow(taux,4) +c3*pow(taux,3) + c2*pow(taux,2) + c1*taux + c0;
          polyy = c5*pow(tauy,5) + c4*pow(tauy,4) +c3*pow(tauy,3) + c2*pow(tauy,2) + c1*tauy + c0;

          polyx = qx<=q1 ? 0. : polyx;  polyx = qx>=q2 ? 1. : polyx;
          polyy = qy<=q1 ? 0. : polyy;  polyy = qy>=q2 ? 1. : polyy;
           
          if (type==100) { //X Pml
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*polyx;
          } else if (type==200) { //Y Pml
             mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*polyy;
          } else if (type==300) { //XY Pml
            if(x>xmax || x<xmin)
             mesh->pmlSigmaX[Nnodes*pmlId + n] = xsigma*polyx;
            if(y>ymax || y<ymin)
             mesh->pmlSigmaY[Nnodes*pmlId + n] = ysigma*polyy;
          }
        }

      }

    }

    printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
    mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

    //mesh->Nfields = 6;

    if(strstr(options,"SAAB")  || strstr(options,"SRAB") ){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqy);
    }



    if(strstr(options,"LSERK")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlresqx = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlresqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqx);
      mesh->o_pmlresqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlresqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqy);
      mesh->o_pmlresqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlresqy);
    }


    if(strstr(options,"SARK")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqx = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));

      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlrhsqy = (dfloat*) calloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qSx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_pmlrhsqx  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qSy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_pmlrhsqy  = mesh->device.malloc(mesh->Nrhs*mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlrhsqy);
    }


    if(strstr(options,"LSIMEX")){
      mesh->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      mesh->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->Nfields, sizeof(dfloat));
      // set up PML on DEVICE 
      mesh->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);   
      mesh->o_qSx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qYx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);
      mesh->o_qZx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqx);

      mesh->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qSy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qYy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
      mesh->o_qZy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->Nfields*sizeof(dfloat), mesh->pmlqy);
    }

    if (mesh->pmlNelements) {
      mesh->o_pmlSigmaX     = mesh->device.malloc(mesh->pmlNelements*Nnodes*sizeof(dfloat),mesh->pmlSigmaX);
      mesh->o_pmlSigmaY     = mesh->device.malloc(mesh->pmlNelements*Nnodes*sizeof(dfloat),mesh->pmlSigmaY);
      mesh->o_pmlElementIds = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlElementIds);
      mesh->o_pmlIds        = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlIds);
    }
      
   free(pmlIds);
 //    //
 //    printf("Wwitting Pml profile\n" );
 //     iint fid = 0; 
 // for (iint es=0;es<mesh->pmlNelements;es++){
 //       iint e     = mesh->pmlElementIds[es];
 //       iint pmlId = mesh->pmlIds[es];   
 //    for(iint n=0;n<Nnodes;n++){
 //      mesh->q[mesh->Nfields*(n + e*mesh->Np) + 0] = mesh->pmlSigmaX[Nnodes*pmlId + n];
 //      mesh->q[mesh->Nfields*(n + e*mesh->Np) + 1] = mesh->pmlSigmaY[Nnodes*pmlId + n];
 //    }
  // }
  #if 0
   dfloat *cubProjectT = (dfloat*) calloc(mesh->cubNp*mesh->Np, sizeof(dfloat));
   mesh->o_cubProjectT.copyTo(cubProjectT);
    // free(pmlIds);
     iint fid = 0; 
 for (iint es=0;es<mesh->pmlNelements;es++){
       iint e     = mesh->pmlElementIds[es];
       iint pmlId = mesh->pmlIds[es];   
    for(iint n=0;n<mesh->Np;n++)
    {
      dfloat q1 = 0; dfloat q2 = 0; 
      for(iint m=0; m<mesh->cubNp;++m){
         dfloat prj = cubProjectT[m*mesh->Np+n];
         q1 += prj*mesh->pmlSigmaX[Nnodes*pmlId + m];
         q2 += prj*mesh->pmlSigmaY[Nnodes*pmlId + m];
      }

      mesh->q[mesh->Nfields*(n + e*mesh->Np) + 0] = q1;
      mesh->q[mesh->Nfields*(n + e*mesh->Np) + 1] = q2;
    }
  }
  

  char fname[BUFSIZ];
  sprintf(fname, "pmlProfile1.vtu");
  boltzmannPlotVTU2D(mesh, fname);

   //
  printf("NonPmlElements: %d  PmlElements: %d \n", mesh->nonPmlNelements, mesh->pmlNelements);

  
  for (iint es=0;es<mesh->nonPmlNelements;es++){
       iint e      = mesh->nonPmlElementIds[es];
       iint id = e*mesh->Nverts;  
        dfloat x1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat x2 = mesh->EX[id+1];
        dfloat x3 = mesh->EX[id+2];

        dfloat y1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat y2 = mesh->EY[id+1];
        dfloat y3 = mesh->EY[id+2];
       printf("Element: %d  x1: %.5f x2: %.5f x3: %.5f y1: %.5f y2: %.5f y3: %.5f \n", e, x1,x2,x3,y1,y2,y3);
  }   


  for (iint es=0;es<mesh->pmlNelements;es++){
       iint e     = mesh->pmlElementIds[es];
       iint pmlId = mesh->pmlIds[es];
       iint id = e*mesh->Nverts;  
       //
       dfloat x1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat x2 = mesh->EX[id+1];
        dfloat x3 = mesh->EX[id+2];

        dfloat y1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat y2 = mesh->EY[id+1];
        dfloat y3 = mesh->EY[id+2];
       printf("Element: %d   PmlId = %d x1: %.5f x2: %.5f x3: %.5f y1: %.5f y2: %.5f y3: %.5f\n", e, pmlId, x1,x2,x3,y1,y2,y3);
  }   



  #endif
  
  }

}
