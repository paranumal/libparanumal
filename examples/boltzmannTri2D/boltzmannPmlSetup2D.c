#include "boltzmann2D.h"

void boltzmannPmlSetup2D(bns_t *bns, setupAide &options){

  mesh2D *mesh = bns->mesh; 

  
  dfloat xsigma    = 100., ysigma  = 100.;
  
  dfloat xsmin    = -2.0, xsmax = 2.0; // map x to this range to control erf profile 
  dfloat ysmin    = -2.0, ysmax = 2.0; // map y to this range to control erf profile 

  dfloat q1       = 0.0, q2 = sqrt(2);    // Ramped profile coefficients 0<q1, 1>q2 , q1<exp<q2 or polynomial
  dfloat alphamax = 0.2; 
  dfloat betamax  = 1.0; 
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

  //printf(" Sigma Scale: %d\n",mesh->Ntscale);
  //count the pml elements
  mesh->pmlNelements=0;


  for (int m=0;m<mesh->Nelements;m++) {
    int e = mesh->nonPmlElementIds[m];
    int type = mesh->elementInfo[e];
    if ((type==100)||(type==200)||(type==300)) {
      mesh->pmlNelements++;
    }
  }

  //set up the pml
  if (mesh->pmlNelements) {

    //construct a numbering of the pml elements
    int *pmlIds = (int *) calloc(mesh->Nelements,sizeof(int));
    int pmlcnt  = 0;
    int nonpmlcnt = 0;
    for (int e=0;e<mesh->Nelements;e++) {
      int type = mesh->elementInfo[e];
      if ((type==100)||(type==200)||(type==300))  //pml element
        pmlIds[e] = pmlcnt++;
    }

    //set up lists of pml elements and remove the pml elements from the nonpml element list

    mesh->pmlElementIds = (int *) calloc(mesh->pmlNelements,sizeof(int));
    mesh->pmlIds        = (int *) calloc(mesh->pmlNelements,sizeof(int));


    pmlcnt = 0;
    nonpmlcnt = 0;
    for (int m=0;m<mesh->Nelements;m++){
      int e = mesh->nonPmlElementIds[m];
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
    mesh->nonPmlElementIds = (int*) realloc(mesh->nonPmlElementIds,mesh->nonPmlNelements*sizeof(int));


    int Nnodes = 0;
    if(options.compareArgs("RELAXATION TYPE","CUBATURE")){ // !!!!!!!!!!!!!!!
      //
      printf("Setting PML Coefficient for Cubature Integration\n");
      //set up damping parameter
      bns->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat));
      bns->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->cubNp,sizeof(dfloat));  
      Nnodes = mesh->cubNp;    
    }
     else{
      printf("Setting PML Coefficients for Nodal Collocation Integration\n");
      //set up damping parameter
      bns->pmlSigmaX = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));
      bns->pmlSigmaY = (dfloat *) calloc(mesh->pmlNelements*mesh->Np,sizeof(dfloat));
      Nnodes = mesh->Np; 
    }
   
    //find the bounding box of the whole domain and interior domain
    dfloat xmin = 1e9, xmax =-1e9;
    dfloat ymin = 1e9, ymax =-1e9;
    dfloat pmlxmin = 1e9, pmlxmax =-1e9;
    dfloat pmlymin = 1e9, pmlymax =-1e9;
    for (int e=0;e<mesh->Nelements;e++) {
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



    int order = 2; 
    //
    if(options.compareArgs("PML PROFILE","CONSTANT"))
      order = 0; 
    if(options.compareArgs("PML PROFILE","LINEAR"))
      order = 1; 
    if(options.compareArgs("PML PROFILE","QUADRATIC"))
      order = 2;
    if(options.compareArgs("PML PROFILE","FORTHORDER"))
      order = 4;
    if(options.compareArgs("PML PROFILE","ERFFUNCTION"))
      order =1;

    dfloat xmaxScale = pow(pmlxmax-xmax,order);
    dfloat xminScale = pow(pmlxmin-xmin,order);
    dfloat ymaxScale = pow(pmlymax-ymax,order);
    dfloat yminScale = pow(pmlymin-ymin,order);


    printf("xmaxScale = %.10e \n", xmaxScale);

    //set up the damping factor
    for (int es=0;es<mesh->pmlNelements;es++){
        int e     = mesh->pmlElementIds[es];
        int pmlId = mesh->pmlIds[es];
        int type  = mesh->elementInfo[e];
        int id    = e*mesh->Nverts;

        dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = mesh->EX[id+1];
        dfloat xe3 = mesh->EX[id+2];

        dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = mesh->EY[id+1];
        dfloat ye3 = mesh->EY[id+2];

        for(int n=0;n<Nnodes;++n){ /* for each node */
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
          
          // if(!strstr(options,"ERFUNCTION") && !strstr(options,"RAMPFUNCTION") && !strstr(options,"SMOOTHPOLYNOMIAL")){
          if (type==100) { //X Pml
            if(x>xmax)
              bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmin,order)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymin,order)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*pow(x-xmin,order)/xminScale;
            if(y>ymax)
              bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*pow(y-ymin,order)/yminScale;
          }
        // }

        // else if(strstr(options,"ERFUNCTION")){
        //   if (type==100) { //X Pml
        //     if(x>xmax)
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmax)/xmaxScale * (xsmax-xsmin)));
        //     if(x<xmin)
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmin)/xminScale * (xsmax-xsmin)));
        //   } else if (type==200) { //Y Pml
        //     if(y>ymax)
        //       bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymax)/ymaxScale * (ysmax-ysmin)));
        //     if(y<ymin)
        //       bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymin)/yminScale * (ysmax-ysmin)));
        //   } else if (type==300) { //XY Pml
        //     if(x>xmax)
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmax)/xmaxScale * (xsmax-xsmin)));
        //     if(x<xmin)
        //       bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*0.5*(1.0+erf(xsmin + (x-xmin)/xminScale * (xsmax-xsmin)));
        //     if(y>ymax)
        //       bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymax)/ymaxScale * (ysmax-ysmin)));
        //     if(y<ymin)
        //       bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*0.5*(1.0+erf(ysmin + (y-ymin)/yminScale * (ysmax-ysmin)));
        //   }
        // }

        // else if(strstr(options,"RAMPFUNCTION")){

        //   dfloat qx=0,      qy = 0;
        //   dfloat taux = 0,  tauy = 0;
        //   dfloat rampx = 0, rampy = 0;
          
        //   // stretch x to [0 1] interval
        //   if (type==100) { //X Pml
        //     if(x>xmax)
        //       qx   = (x-xmax)/(pmlxmax-xmax);
        //     if(x<xmin)
        //       qx    = (x-xmin)/(pmlxmin-xmin);
        //   } else if (type==200) { //Y Pml
        //     if(y>ymax)
        //       qy    = (y-ymax)/(pmlymax-ymax);
        //     if(y<ymin)
        //       qy    = (y-ymin)/(pmlymin-ymin);
        //   } else if (type==300) { //XY Pml
        //     if(x>xmax)
        //      qx   = (x-xmax)/(pmlxmax-xmax);
        //     if(x<xmin)
        //      qx    = (x-xmin)/(pmlxmin-xmin);
        //     if(y>ymax)
        //      qy    = (y-ymax)/(pmlymax-ymax);
        //     if(y<ymin)
        //      qy    = (y-ymin)/(pmlymin-ymin);
        //   }
          
        //   // Renormalize for ramp profile i.e. tau = [0,1]
        //   taux  =  (qx-q1)/(q2-q1);
        //   tauy  =  (qy-q1)/(q2-q1);

        //   rampx = 1.0 -exp(2.*pow(exp(1.0),(-1. /taux))/(taux-1.0));
        //   rampy = 1.0 -exp(2.*pow(exp(1.0),(-1. /tauy))/(tauy-1.0));
        //   rampx = qx<=q1 ? 0. : rampx;  rampx = qx>=q2 ? 1. : rampx;
        //   rampy = qy<=q1 ? 0. : rampy;  rampy = qy>=q2 ? 1. : rampy;
           
        //   if (type==100) { //X Pml
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*rampx;
        //   } else if (type==200) { //Y Pml
        //      bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*rampy;
        //   } else if (type==300) { //XY Pml
        //     if(x>xmax || x<xmin)
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*rampx;
        //     if(y>ymax || y<ymin)
        //      bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*rampy;
        //   }
        // }

        // else if(strstr(options,"SMOOTHPOLYNOMIAL")){

        //   dfloat qx=0,      qy = 0;
        //   dfloat taux = 0,  tauy = 0;
        //   dfloat polyx = 0, polyy = 0;

        //   if (type==100) { //X Pml
        //     if(x>xmax)
        //       qx   = (x-xmax)/(pmlxmax-xmax);
        //     if(x<xmin)
        //       qx    = (x-xmin)/(pmlxmin-xmin);
        //   } else if (type==200) { //Y Pml
        //     if(y>ymax)
        //       qy    = (y-ymax)/(pmlymax-ymax);
        //     if(y<ymin)
        //       qy    = (y-ymin)/(pmlymin-ymin);
        //   } else if (type==300) { //XY Pml
        //     if(x>xmax)
        //      qx   = (x-xmax)/(pmlxmax-xmax);
        //     if(x<xmin)
        //      qx    = (x-xmin)/(pmlxmin-xmin);
        //     if(y>ymax)
        //      qy    = (y-ymax)/(pmlymax-ymax);
        //     if(y<ymin)
        //      qy    = (y-ymin)/(pmlymin-ymin);
        //   }
        //   // Third Order
        //   taux  =  (qx-q1)/(q2-q1);
        //   tauy  =  (qy-q1)/(q2-q1);
        //   // polyx = c3*pow(taux,3) + c2*pow(taux,2) +c1*taux +c0;
        //   // polyy = c3*pow(tauy,3) + c2*pow(tauy,2) +c1*tauy +c0;
        //   // // Fifth Order
        //   polyx = c5*pow(taux,5) + c4*pow(taux,4) +c3*pow(taux,3) + c2*pow(taux,2) + c1*taux + c0;
        //   polyy = c5*pow(tauy,5) + c4*pow(tauy,4) +c3*pow(tauy,3) + c2*pow(tauy,2) + c1*tauy + c0;

        //   polyx = qx<=q1 ? 0. : polyx;  polyx = qx>=q2 ? 1. : polyx;
        //   polyy = qy<=q1 ? 0. : polyy;  polyy = qy>=q2 ? 1. : polyy;
           
        //   if (type==100) { //X Pml
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*polyx;
        //   } else if (type==200) { //Y Pml
        //      bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*polyy;
        //   } else if (type==300) { //XY Pml
        //     if(x>xmax || x<xmin)
        //      bns->pmlSigmaX[Nnodes*pmlId + n] = xsigma*polyx;
        //     if(y>ymax || y<ymin)
        //      bns->pmlSigmaY[Nnodes*pmlId + n] = ysigma*polyy;
        //   }


        // }

      }

    }


 #if 0   
    // Setup Pml Parameter beta and initiliaze to zero 
    mesh->pmlBetaX = (dfloat *) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
    mesh->pmlBetaY = (dfloat *) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

    dfloat betap =   1.0; 
    dfloat betam =  -1.0; // was -1.0

    dfloat betaxMax = sqrt(3.0);  // was 1.0
    dfloat betayMax = sqrt(3.0); // was 1.0
   
    for (int es=0;es<mesh->pmlNelements;es++){
        int e     = mesh->pmlElementIds[es];
        int pmlId = mesh->pmlIds[es];
        int type  = mesh->elementInfo[e];
        int id    = e*mesh->Nverts;

        dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
        dfloat xe2 = mesh->EX[id+1];
        dfloat xe3 = mesh->EX[id+2];

        dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
        dfloat ye2 = mesh->EY[id+1];
        dfloat ye3 = mesh->EY[id+2];

        for(int n=0;n<mesh->Np;++n){ /* for each node */
            dfloat x = mesh->x[n + e*mesh->Np];
            dfloat y = mesh->y[n + e*mesh->Np];
            //
          if (type==100) { //X Pml
            if(x>xmax)
              mesh->pmlBetaX[mesh->Np*e + n] = betap*betaxMax*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              mesh->pmlBetaX[mesh->Np*e + n] = betam*betaxMax*pow(x-xmin,order)/xminScale;
          } else if (type==200) { //Y Pml
            if(y>ymax)
              mesh->pmlBetaY[mesh->Np*e + n] = betap*betayMax*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              mesh->pmlBetaY[mesh->Np*e + n] = betam*betayMax*pow(y-ymin,order)/yminScale;
          } else if (type==300) { //XY Pml
            if(x>xmax)
              mesh->pmlBetaX[mesh->Np*e + n] = betap*betaxMax*pow(x-xmax,order)/xmaxScale;
            if(x<xmin)
              mesh->pmlBetaX[mesh->Np*e + n] = betam*betaxMax*pow(x-xmin,order)/xminScale;
            if(y>ymax)
              mesh->pmlBetaY[mesh->Np*e + n] = betap*betayMax*pow(y-ymax,order)/ymaxScale;
            if(y<ymin)
              mesh->pmlBetaY[mesh->Np*e + n] = betam*betayMax*pow(y-ymin,order)/yminScale;
          }  


          

        }
      }

      mesh->o_pmlBetaX      = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),mesh->pmlBetaX);
      mesh->o_pmlBetaY      = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat),mesh->pmlBetaY); 

#endif

    printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
    mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

    //bns->Nfields = 6;

    if(options.compareArgs("TIME INTEGRATOR","SRSAAB")  || options.compareArgs("TIME INTEGRATOR","SRAB") ){
      bns->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx = (dfloat*) calloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy = (dfloat*) calloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
    }



    if(options.compareArgs("TIME INTEGRATOR","LSERK")){
      bns->pmlqx     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlresqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqy     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlresqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);
      bns->o_pmlresqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlresqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
      bns->o_pmlresqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlresqy);
    }

    if(options.compareArgs("TIME INTEGRATOR","DOPRI5") || options.compareArgs("TIME INTEGRATOR","XDOPRI") || 
       options.compareArgs("TIME INTEGRATOR","SAADRK") ){
      bns->pmlqx     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      
      bns->pmlqy     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkqx      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkqy      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkrhsqx   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkrhsqy   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    
      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
      
      bns->o_rkqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqx);
      bns->o_rkqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqy);
      
      bns->o_rkrhsqx  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqx);
      bns->o_rkrhsqy  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqy);

    }


    if(options.compareArgs("TIME INTEGRATOR","SARK")){
      bns->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx = (dfloat*) calloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy = (dfloat*) calloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));


      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_qSx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_qSy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(bns->Nrhs*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
    }


    if(options.compareArgs("TIME INTEGRATOR","LSIMEX")){
      bns->pmlqx    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlqy    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
 
      // set up PML on DEVICE 
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);   
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);
      bns->o_qYx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      // mesh->o_qZx       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
      bns->o_qYy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      // mesh->o_qZy       = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
    }

    if(options.compareArgs("TIME INTEGRATOR","IMEXRK")){

      // printf("Here initializing IMEXRK PML\n");

      bns->pmlqx     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqx  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      
      bns->pmlqy     = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->pmlrhsqy  = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkqx      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkqy      = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));

      bns->rkrhsqx   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
      bns->rkrhsqy   = (dfloat*) calloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields, sizeof(dfloat));
    
      // set up PML on DEVICE    
      bns->o_pmlqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqx);
      bns->o_pmlrhsqx  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqx);

      bns->o_pmlqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlqy);
      bns->o_pmlrhsqy  = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->pmlrhsqy);
      
      bns->o_rkqx     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqx);
      bns->o_rkqy     = mesh->device.malloc(mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkqy);
      
      bns->o_rkrhsqx  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqx);
      bns->o_rkrhsqy  = mesh->device.malloc(bns->NrkStages*mesh->pmlNelements*mesh->Np*bns->Nfields*sizeof(dfloat), bns->rkrhsqy);
    }

    if (mesh->pmlNelements) {
      bns->o_pmlSigmaX     = mesh->device.malloc(mesh->pmlNelements*Nnodes*sizeof(dfloat),bns->pmlSigmaX);
      bns->o_pmlSigmaY     = mesh->device.malloc(mesh->pmlNelements*Nnodes*sizeof(dfloat),bns->pmlSigmaY);
      mesh->o_pmlElementIds = mesh->device.malloc(mesh->pmlNelements*sizeof(int), mesh->pmlElementIds);
      mesh->o_pmlIds        = mesh->device.malloc(mesh->pmlNelements*sizeof(int), mesh->pmlIds);
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
 //    for(int n=0;n<Nnodes;n++){
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
 //         q1 += prj*bns->pmlSigmaX[Nnodes*pmlId + m];
 //         q2 += prj*bns->pmlSigmaY[Nnodes*pmlId + m];
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
