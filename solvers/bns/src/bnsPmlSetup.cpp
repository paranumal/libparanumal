/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "bns.hpp"

void bns_t::PmlSetup(){

  //build pml element lists
  mesh.PmlSetup();

  pmlcubature=0;

  //set up the pml
  if (mesh.NpmlElements) {
    //get settings from solver
    settings.getSetting("PML PROFILE ORDER", pmlOrder);
    settings.getSetting("PML SIGMAX MAX", sigmaXmax);
    settings.getSetting("PML SIGMAY MAX", sigmaYmax);
    settings.getSetting("PML SIGMAZ MAX", sigmaZmax);
    pmlcubature = (settings.compareSetting("PML INTEGRATION", "CUBATURE")) ? 1:0;

    pmlAlpha = 0.2; //hardcoded

    //find the bounding box of the whole domain and interior domain
    dfloat xmin = 1e9, xmax =-1e9;
    dfloat ymin = 1e9, ymax =-1e9;
    dfloat zmin = 1e9, zmax =-1e9;
    dfloat pmlxmin = 1e9, pmlxmax =-1e9;
    dfloat pmlymin = 1e9, pmlymax =-1e9;
    dfloat pmlzmin = 1e9, pmlzmax =-1e9;

    //for all elements
    for (dlong e=0;e<mesh.Nelements;e++) {
      for (int n=0;n<mesh.Nverts;n++) {
        dfloat x =0, y=0, z = 0;
        x  = mesh.EX[e*mesh.Nverts+n];
        y  = mesh.EY[e*mesh.Nverts+n];

        pmlxmin = (pmlxmin > x) ? x : pmlxmin;
        pmlymin = (pmlymin > y) ? y : pmlymin;
        pmlxmax = (pmlxmax < x) ? x : pmlxmax;
        pmlymax = (pmlymax < y) ? y : pmlymax;

        if(mesh.dim==3){
          z = mesh.EZ[e*mesh.Nverts+n];
          pmlzmin = (pmlzmin > z) ? z : pmlzmin;
          pmlzmax = (pmlzmax < z) ? z : pmlzmax;
        }
      }
    }

    //for interior elements
    for (dlong m=0;m<mesh.NnonPmlElements;m++) {
      dlong e = mesh.nonPmlElements[m];
      for (int n=0;n<mesh.Nverts;n++) {
        dfloat x = 0., y = 0., z = 0.;
        x = mesh.EX[e*mesh.Nverts+n];
        y = mesh.EY[e*mesh.Nverts+n];
        xmin = (xmin > x) ? x : xmin;
        ymin = (ymin > y) ? y : ymin;
        xmax = (xmax < x) ? x : xmax;
        ymax = (ymax < y) ? y : ymax;
        if(mesh.dim==3){
          z = mesh.EZ[e*mesh.Nverts+n];
          zmin = (zmin > z) ? z : zmin;
          zmax = (zmax < z) ? z : zmax;
        }
      }
    }

    dfloat xmaxScale =0., xminScale=0.;
    dfloat ymaxScale =0., yminScale=0.;
    dfloat zmaxScale =0., zminScale=0.;

    xmaxScale = pow(pmlxmax-xmax,pmlOrder);
    xminScale = pow(pmlxmin-xmin,pmlOrder);
    ymaxScale = pow(pmlymax-ymax,pmlOrder);
    yminScale = pow(pmlymin-ymin,pmlOrder);
    if(mesh.dim==3){
      zmaxScale = pow(pmlzmax-zmax,pmlOrder);
      zminScale = pow(pmlzmin-zmin,pmlOrder);
    }

    // Set the size of pml nodes
    int pmlNp = (pmlcubature) ? mesh.cubNp : mesh.Np;
    int pmlNq = (pmlcubature) ? mesh.cubNq : mesh.Nq;

    dfloat *pmlr, *pmls, *pmlt;
    if(pmlcubature){
      pmlr = mesh.cubr;
      pmls = mesh.cubs;
      pmlt = mesh.cubt;
    } else {
      pmlr = mesh.r;
      pmls = mesh.s;
      pmlt = mesh.t;
    }

    // printf("Setting PML Coefficient \n");
    //set up damping parameter
    pmlSigma = (dfloat *) calloc(mesh.dim*mesh.NpmlElements*pmlNp,sizeof(dfloat));

    for (dlong m=0;m<mesh.NpmlElements;m++){
      dlong e     = mesh.pmlElements[m];
      hlong type  = mesh.elementInfo[e];

      //element vertices
      const dfloat *xe = mesh.EX + e*mesh.Nverts;
      const dfloat *ye = mesh.EY + e*mesh.Nverts;
      const dfloat *ze = mesh.EZ + e*mesh.Nverts;

      for(int n=0;n<pmlNp;++n){ /* for each node */
        dfloat x  = 0, y  = 0, z  = 0;
        dfloat rn = 0, sn = 0, tn = 0;
        if(mesh.elementType==TRIANGLES){
          rn = pmlr[n];
          sn = pmls[n];

          x = -0.5*(rn+sn)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2];
          y = -0.5*(rn+sn)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2];
        } else if(mesh.elementType==QUADRILATERALS){
          const int i = n%pmlNq;
          const int j = n/pmlNq;
          rn = pmlr[i];
          sn = pmlr[j];

          x =  0.25*( (1.0-rn)*(1-sn)*xe[0]+(1.0-rn)*(1+sn)*xe[1]+(1.0+rn)*(1+sn)*xe[2]+(1.0+rn)*(1-sn)*xe[3]);
          y =  0.25*( (1.0-rn)*(1-sn)*ye[0]+(1.0-rn)*(1+sn)*ye[1]+(1.0+rn)*(1+sn)*ye[2]+(1.0+rn)*(1-sn)*ye[3]);
        } else if(mesh.elementType==TETRAHEDRA){
          rn = pmlr[n];
          sn = pmls[n];
          tn = pmlt[n];

          x = -0.5*(rn+sn+tn+1)*xe[0] + 0.5*(1+rn)*xe[1] + 0.5*(1+sn)*xe[2] + 0.5*(tn+1)*xe[3];
          y = -0.5*(rn+sn+tn+1)*ye[0] + 0.5*(1+rn)*ye[1] + 0.5*(1+sn)*ye[2] + 0.5*(tn+1)*ye[3];
          z = -0.5*(rn+sn+tn+1)*ze[0] + 0.5*(1+rn)*ze[1] + 0.5*(1+sn)*ze[2] + 0.5*(tn+1)*ze[3];
        } else if(mesh.elementType==HEXAHEDRA){
          const int i = n%pmlNq;
          const int j = (n/pmlNq)%pmlNq;
          const int k = (n/pmlNq)/pmlNq;
          rn = pmlr[i];
          sn = pmlr[j];
          tn = pmlr[k];

          x =
            +0.125*(1-rn)*(1-sn)*(1-tn)*xe[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*xe[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*xe[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*xe[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*xe[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*xe[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*xe[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*xe[7];

          y =
            +0.125*(1-rn)*(1-sn)*(1-tn)*ye[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*ye[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*ye[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*ye[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*ye[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*ye[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*ye[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*ye[7];

          z =
            +0.125*(1-rn)*(1-sn)*(1-tn)*ze[0]
            +0.125*(1+rn)*(1-sn)*(1-tn)*ze[1]
            +0.125*(1+rn)*(1+sn)*(1-tn)*ze[2]
            +0.125*(1-rn)*(1+sn)*(1-tn)*ze[3]
            +0.125*(1-rn)*(1-sn)*(1+tn)*ze[4]
            +0.125*(1+rn)*(1-sn)*(1+tn)*ze[5]
            +0.125*(1+rn)*(1+sn)*(1+tn)*ze[6]
            +0.125*(1-rn)*(1+sn)*(1+tn)*ze[7];
        }

        if (type==100) { //X Pml
          if(x>xmax) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmax,pmlOrder)/xmaxScale;
          if(x<xmin) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmin,pmlOrder)/xminScale;
        } else if (type==200) { //Y Pml
          if(y>ymax) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymax,pmlOrder)/ymaxScale;
          if(y<ymin) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymin,pmlOrder)/yminScale;
        } else if (type==300) { //XY Pml
          if(x>xmax) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmax,pmlOrder)/xmaxScale;
          if(x<xmin) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmin,pmlOrder)/xminScale;
          if(y>ymax) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymax,pmlOrder)/ymaxScale;
          if(y<ymin) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymin,pmlOrder)/yminScale;
        }

        if(mesh.dim==3){
          if (type==400) { //Z Pml
            if(z>zmax) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmax,pmlOrder)/zmaxScale;
            if(z<zmin) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmin,pmlOrder)/zminScale;
          } else if (type==500) {//XZ Pml
            if(x>xmax) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmax,pmlOrder)/xmaxScale;
            if(x<xmin) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmin,pmlOrder)/xminScale;
            if(z>zmax) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmax,pmlOrder)/zmaxScale;
            if(z<zmin) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmin,pmlOrder)/zminScale;
          } else if (type==600){ //YZ Pml
            if(y>ymax) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymax,pmlOrder)/ymaxScale;
            if(y<ymin) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymin,pmlOrder)/yminScale;
            if(z>zmax) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmax,pmlOrder)/zmaxScale;
            if(z<zmin) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmin,pmlOrder)/zminScale;
          } else if (type==700){ //XYZ Pml
            if(x>xmax) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmax,pmlOrder)/xmaxScale;
            if(x<xmin) pmlSigma[mesh.dim*pmlNp*m + 0*pmlNp + n] = sigmaXmax*pow(x-xmin,pmlOrder)/xminScale;
            if(y>ymax) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymax,pmlOrder)/ymaxScale;
            if(y<ymin) pmlSigma[mesh.dim*pmlNp*m + 1*pmlNp + n] = sigmaYmax*pow(y-ymin,pmlOrder)/yminScale;
            if(z>zmax) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmax,pmlOrder)/zmaxScale;
            if(z<zmin) pmlSigma[mesh.dim*pmlNp*m + 2*pmlNp + n] = sigmaZmax*pow(z-zmin,pmlOrder)/zminScale;
          }
        }
      }
    }

    // printf("# of PML elements: %d and # of Non-PML elements: %d \n",mesh.NpmlElements, mesh.Nelements-mesh.NpmlElements);
    if (mesh.NpmlElements)
      o_pmlSigma = mesh.device.malloc(mesh.dim*mesh.NpmlElements*pmlNp*sizeof(dfloat),pmlSigma);
  }
}
