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

#include "esdg.hpp"

void esdg_t::SetupTri2D(){

  settings.getSetting("FLUX DEGREE INCREMENT", fluxN);
  fluxN += mesh.N;
  fluxNp = (mesh.dim==2) ? (fluxN+1)*(fluxN+2)/2 : (fluxN+1)*(fluxN+2)*(fluxN+3)/6; 

    // use this for outflow viscosity ramp
  dfloat outflowMu, outflowXmin, outflowXmax, meanMu;
  settings.getSetting("LAME MU", meanMu);
  settings.getSetting("OUTFLOW MU", outflowMu);
  
  settings.getSetting("OUTFLOW LAYER XMIN", outflowXmin);
  settings.getSetting("OUTFLOW LAYER XMAX", outflowXmax);



  
  memory<dfloat> fluxr, fluxs, fluxt, fluxVflux, fluxVsoln, fluxMMflux, fluxIsoln, fluxIsolnT, fluxFfluxT, fluxMeanWeights;

  fluxr.malloc(fluxNp);
  fluxs.malloc(fluxNp);
  fluxt.malloc(fluxNp);

  fluxVflux .malloc(fluxNp*fluxNp);
  fluxVsoln .malloc(mesh.Np*fluxNp);

  fluxMMflux.malloc(fluxNp*fluxNp);
  fluxIsoln .malloc(fluxNp*mesh.Np);

  fluxIsolnT.malloc(fluxNp*mesh.Np);
  fluxFfluxT.malloc(fluxNp*fluxNp);

  fluxMeanWeights.malloc(fluxNp);
  
  
  //    mesh.ChebyshevNodesTri2D(fluxN, fluxr, fluxs);
  mesh.NodesTri2D(fluxN, fluxr, fluxs);
  
  mesh.VandermondeTri2D(fluxN,  fluxr,    fluxs,    fluxVflux);
  mesh.VandermondeTri2D(mesh.N,       fluxr,    fluxs,    fluxVsoln);
  mesh.MassMatrixTri2D(fluxNp, fluxVflux, fluxMMflux);
  
  // interpolates from solution nodes to flux
  mesh.InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, fluxr, fluxs, fluxIsoln);

  // compute weights used for computing mean  on flux nodes
  for(int n=0;n<fluxNp;++n){
    for(int m=0;m<fluxNp;++m){
      fluxMeanWeights[n] += fluxMMflux[n*fluxNp+m];
    }
  }

  // store soln to flux interpolation matrix
  for(int n=0;n<fluxNp;++n){
    for(int m=0;m<mesh.Np;++m){
      fluxIsolnT[n + m*fluxNp] = fluxIsoln[n*mesh.Np+m];
    }
  }
  
  // built project from flux nodes to super nodes

  // filter from flux to flux nodes thorough solution space
  for(int n=0;n<fluxNp;++n){
    for(int m=0;m<fluxNp;++m){
      dfloat resFilt = 0;
      for(int j=0;j<mesh.Np;++j){
	for(int i=0;i<fluxNp;++i){
	  resFilt += fluxVsoln[mesh.Np*n + j]*fluxVsoln[mesh.Np*i+j]*fluxMMflux[i*fluxNp+m];
	}
      }
      fluxFfluxT[n+m*fluxNp] = resFilt;
    }
  }

  o_fluxIsolnT = platform.malloc<dfloat>(fluxNp*mesh.Np, fluxIsolnT);
  o_fluxFfluxT = platform.malloc<dfloat>(fluxNp*fluxNp, fluxFfluxT);
  o_fluxMMfluxT = platform.malloc<dfloat>(fluxNp*fluxNp, fluxMMflux);
  o_fluxMeanWeights = platform.malloc<dfloat>(fluxNp, fluxMeanWeights);

  // set up cubature grid
  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  entropyStable  = (settings.compareSetting("ADVECTION TYPE", "ENTROPYSTABLE")) ? 1:0;

  mesh.CubatureSetup();
  mesh.CubaturePhysicalNodes();

  // TURNED OFF FOR A MOMENT
  Nfields   = mesh.dim + 2;
  Ngrads = mesh.dim*mesh.dim; // fix this
  
  // WATCH OUT - THIS DISABLES VARIABLE FILTERING
  reconN = mesh.N; // -reconN; // WARNING
  reconNp = (mesh.dim==2) ? ((reconN+1)*(reconN+2))/2 :
    ((reconN+1)*(reconN+2)*(reconN+3)/6);

  // build element centers
  memory<dfloat> cx, cy, cz;
  cx.malloc(mesh.Nelements);
  cy.malloc(mesh.Nelements);
  cz.malloc(mesh.Nelements);

  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat cxe = 0, cye = 0, cze = 0;
    for(int n=0;n<mesh.Nverts;++n){
      cxe += mesh.EX[e*mesh.Nverts+n];
      cye += mesh.EY[e*mesh.Nverts+n];
      if(mesh.dim==3)
	cze += mesh.EZ[e*mesh.Nverts+n];
    }
    cx[e] = cxe/mesh.Nverts;
    cy[e] = cye/mesh.Nverts;
    cz[e] = cze/mesh.Nverts;
  }
  o_cx = platform.malloc<dfloat>(mesh.Nelements, cx);
  o_cy = platform.malloc<dfloat>(mesh.Nelements, cy);
  o_cz = platform.malloc<dfloat>(mesh.Nelements, cz);


  // face cubature
  // build interpolation matrix: nodes to surface nodes
  // HACKITY-HACK (not saved anywhere)
  memory<dfloat> ir, is, it, iw;
  ir.malloc(mesh.intNfp*mesh.Nfaces);
  is.malloc(mesh.intNfp*mesh.Nfaces);
  it.malloc(mesh.intNfp*mesh.Nfaces);
  iw.malloc(mesh.intNfp*mesh.Nfaces);

  memory <dfloat> sInterp, sInterpT;
  sInterp.malloc(mesh.intNfp*mesh.Nfaces*mesh.Np);
  sInterpT.malloc(mesh.intNfp*mesh.Nfaces*mesh.Np);
  
  if(mesh.dim==2){
    for(int n=0;n<mesh.intNfp;++n){
      ir[0*mesh.intNfp + n] =  mesh.intr[n];
      ir[1*mesh.intNfp + n] = -mesh.intr[n];
      ir[2*mesh.intNfp + n] = -1.0;
      
      is[0*mesh.intNfp + n] = -1.0;
      is[1*mesh.intNfp + n] =  mesh.intr[n];
      is[2*mesh.intNfp + n] = -mesh.intr[n];
      
      iw[0*mesh.intNfp + n] =  mesh.intw[n];
      iw[1*mesh.intNfp + n] =  mesh.intw[n];
      iw[2*mesh.intNfp + n] =  mesh.intw[n];

    }

    mesh.InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, ir, is, sInterp);
  }
  
  for(int n=0;n<mesh.intNfp*mesh.Nfaces;++n){
    for(int m=0;m<mesh.Np;++m){
      sInterpT[n+m*mesh.intNfp*mesh.Nfaces] = sInterp[n*mesh.Np+m];
    }
  }
  
  o_sInterpT =
    platform.malloc<dfloat>(mesh.Np*mesh.intNfp*mesh.Nfaces, sInterpT);
  
  o_sQ =
    platform.malloc<dfloat>(mesh.Nelements*mesh.intNfp*mesh.Nfaces*Nfields);

  // ------------------------------------------------------------------------------------------------------------------------------------------
  // entropy stable infrastructure following On discretely entropy conservative and entropy stable discontinuous Galerkin methods, Jesse Chan
  
  int esNq = mesh.cubNp;

  memory<dfloat> esRq, esSq, esWq;
  esRq.malloc(esNq);
  esSq.malloc(esNq);
  esWq.malloc(esNq);
  
  // dim==2
  for(int n=0;n<esNq;++n){
    esRq[n] = mesh.cubr[n];
    esSq[n] = mesh.cubs[n];
    esWq[n] = mesh.cubw[n];
  }

  // esR, esS, esW - combined volume quadrature, surface quadrature nodes and weights
  esNp = esNq + mesh.intNfp*mesh.Nfaces;

  
  printf("ENTROPY STABLE:  N      = %d\n", mesh.N);
  printf("ENTROPY STABLE:  intNfp = %d\n", mesh.intNfp);
  printf("ENTROPY STABLE:  cubN   = %d\n", mesh.cubN); 
  printf("ENTROPY STABLE:  cubNp  = %d\n", mesh.cubNp);
  printf("ENTROPY STABLE:  esNp   = %d\n", esNp);

  int esNf = mesh.Nfaces*mesh.intNfp;

  memory<dfloat> esRf, esSf, esWf;
  esRf.malloc(esNf);
  esSf.malloc(esNf);
  esWf.malloc(esNf);
  
  for(int n=0;n<esNf;++n){
    esRf[n] = ir[n];
    esSf[n] = is[n];
    esWf[n] = iw[n];
  }

  memory<dfloat> esR, esS, esW;
  esR.malloc(esNp);
  esS.malloc(esNp);
  esW.malloc(esNp);

  
  for(int n=0;n<esNq;++n){
    esR[n] = esRq[n];
    esS[n] = esSq[n];
    esW[n] = esWq[n];
  }
  for(int n=0;n<esNf;++n){
    esR[n+esNq] = esRf[n];
    esS[n+esNq] = esSf[n];
    esW[n+esNq] = esWf[n];
  }
  
  // interpolation matrices
  memory<dfloat> esIq, esIqT, esIf;
  esIq.malloc(mesh.Np*esNp);
  esIqT.malloc(mesh.Np*esNp);
  esIf.malloc(mesh.Np*esNf);
      
  
  // Ip = identity

  // modified to interpolate to whole q+f quadrature
  // Iq= Vandermonde2D(N,rp,sp)/V;
  mesh.InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, esRq, esSq, esIq);

  for(int n=0;n<esNq;++n){
    for(int m=0;m<mesh.Np;++m){
      esIqT[n+m*esNq] = esIq[n*mesh.Np+m];
    }
  }
  
  // If = Vandermonde2D(N,rf,sf)/V;
  mesh.InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, esRf, esSf, esIf);

  // build combined interpolator
  memory<dfloat> esIqfT;
  esIqfT.malloc(mesh.Np*esNp);
  
  for(int n=0;n<esNq;++n){
    for(int m=0;m<mesh.Np;++m){
      esIqfT[n+m*esNp] = esIqT[n+m*esNq];
    }
  }
  for(int n=0;n<esNf;++n){
    for(int m=0;m<mesh.Np;++m){
      esIqfT[n+esNq+m*esNp] = esIf[n*mesh.Np+m];
    }
  }
  
  // mass and projection matrices
  //  M = inv(V')*inv(V)
  memory<dfloat> esMM, esInvMM;
  esMM.malloc(mesh.Np*esNq); // hmm
  esInvMM.malloc(mesh.Np*esNq); // hmm
  
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<mesh.Np;++m){
      dfloat resMM = 0, resInvMM = 0;

      for(int j=0;j<esNq;++j){
	// rV*rV'*MM
	resMM += esIq[mesh.Np*j + n]*esWq[j]*esIq[mesh.Np*j+m];
      }
      esMM[n*mesh.Np + m] = resMM;
      esInvMM[n*mesh.Np + m] = resMM;      
    }
  }

  linAlg_t::matrixInverse(mesh.Np, esInvMM);

  // Pq = M\(Vq'*diag(wq)); % J's cancel out
  //    = V*V'*(Vq'*diag(wq));

  memory<dfloat> esPq,esPqT;
  esPq.malloc(mesh.Np*esNq);
  esPqT.malloc(mesh.Np*esNq);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      dfloat resPq = 0;

      for(int j=0;j<mesh.Np;++j){
	resPq += esInvMM[n*mesh.Np+j]*(esIq[m*mesh.Np+j]*esWq[m]);
      }

      esPq[n*esNq + m] = resPq;
      esPqT[n+mesh.Np*m] = resPq;
    }
  }

  memory<dfloat>esItMT;
  esItMT.malloc(mesh.Np*esNq);
  
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      dfloat resItM = 0;

      esItMT[m+esNq*n] = esWq[m]*esIq[m*mesh.Np+n];
    }
  }


  
  // Lf = M\(Vf'*diag(wf));
  memory<dfloat>esLf, esLfT;
  esLf.malloc(mesh.Np*esNf);
  esLfT.malloc(mesh.Np*esNf);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNf;++m){
      dfloat resLf = 0;
      for(int j=0;j<mesh.Np;++j){
	resLf += esInvMM[n*mesh.Np+j]*(esIf[m*mesh.Np+j]*esWf[m]);
      }
      esLf[n*esNf + m] = resLf;
      esLfT[n+m*mesh.Np] = resLf;
    }
  }

  memory<dfloat>esIqfLfT;
  esIqfLfT.malloc(esNp*esNf);

  //  printf("esIqfLfT=[\n");
  for(int n=0;n<esNp;++n){
    for(int m=0;m<esNf;++m){
      dfloat resIqfLf = 0;
      for(int j=0;j<mesh.Np;++j){
	resIqfLf += esIqfT[n+j*esNp]*esLf[j*esNf+m];
      }
      esIqfLfT[n+m*esNp] = resIqfLf;
      //      printf("% g ", resIqfLf);
    }
    //    printf("\n");
  }

  o_esIqfLfT = platform.malloc<dfloat>(esNp*esNf, esIqfLfT);

  memory<dfloat>esIqfDrPqT, esIqfDsPqT;
  esIqfDrPqT.malloc(esNp*esNq);
  esIqfDsPqT.malloc(esNp*esNq);

  
  //  printf("esIqfDr,sPqT=[\n");
  for(int n=0;n<esNp;++n){
    for(int m=0;m<esNq;++m){
      dfloat resIqfDrPq = 0;
      dfloat resIqfDsPq = 0;
      for(int j=0;j<mesh.Np;++j){
	for(int k=0;k<mesh.Np;++k){
	  resIqfDrPq += esIqfT[n+j*esNp]*mesh.Dr[j*mesh.Np+k]*esPq[k*esNq+m];
	  resIqfDsPq += esIqfT[n+j*esNp]*mesh.Ds[j*mesh.Np+k]*esPq[k*esNq+m];
	}
      }
      esIqfDrPqT[n+m*esNp] = resIqfDrPq;
      esIqfDsPqT[n+m*esNp] = resIqfDsPq;
      //      printf("(% g, % g) ", resIqfDrPq, resIqfDsPq);
    }
    //    printf("\n");
  }

  o_esIqfDrPqT = platform.malloc<dfloat>(esNp*esNq, esIqfDrPqT);
  o_esIqfDsPqT = platform.malloc<dfloat>(esNp*esNq, esIqfDsPqT);
  
  // scaled reference normals
  // z = zeros(size(rq1D));
  // e = ones(size(rq1D));
  // nrJ = [-z; e; -e];
  // nsJ = [-e; e; -z];

  memory<dfloat> esNrJ, esNsJ;
  esNrJ.malloc(esNf);
  esNsJ.malloc(esNf);
  for(int n=0;n<esNf;++n){
    int face = n/mesh.intNfp;
    dfloat scr = 0, scs = 0;
    
    if(face==0) { scr =  0.; scs = -1; }
    if(face==1) { scr =  1.; scs =  1; }
    if(face==2) { scr = -1.; scs = 0; }
    
    esNrJ[n] = scr;
    esNsJ[n] = scs;
  }
  
  // "modal" differentiation matrices
  //  Qr_modal = M*Dr;
  //  Qs_modal = M*Ds;
  memory<dfloat> esQrModal, esQsModal;
  esQrModal.malloc(mesh.Np*mesh.Np);
  esQsModal.malloc(mesh.Np*mesh.Np);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<mesh.Np;++m){
      dfloat resQrModal = 0, resQsModal = 0;
      for(int j=0;j<mesh.Np;++j){
	resQrModal += esMM[n*mesh.Np+j]*mesh.Dr[j*mesh.Np+m];
	resQsModal += esMM[n*mesh.Np+j]*mesh.Ds[j*mesh.Np+m];
      }
      esQrModal[n*mesh.Np + m] = resQrModal;
      esQsModal[n*mesh.Np + m] = resQsModal;
    }
  }    
  
  // Qr = Pq'*Qr_modal*Pq;
  // Qs = Pq'*Qs_modal*Pq;
  memory<dfloat> esQr, esQs;
  esQr.malloc(esNq*esNq);
  esQs.malloc(esNq*esNq);

  for(int n=0;n<esNq;++n){
    for(int m=0;m<esNq;++m){
      dfloat resQr = 0, resQs = 0;
      for(int j=0;j<mesh.Np;++j){
	for(int k=0;k<mesh.Np;++k){
	  resQr += esPq[j*esNq+n]*esQrModal[j*mesh.Np+k]*esPq[k*esNq+m];
	  resQs += esPq[j*esNq+n]*esQsModal[j*mesh.Np+k]*esPq[k*esNq+m];
	}
      }
      esQr[n*esNq + m] = resQr;
      esQs[n*esNq + m] = resQs;
    }
  }

  // E = If*Pq;
  memory<dfloat> esE;
  esE.malloc(esNf*esNq);
  for(int n=0;n<esNf;++n){    
    for(int m=0;m<esNq;++m){
      dfloat resE = 0;
      for(int j=0;j<mesh.Np;++j){
	resE += esIf[n*mesh.Np+j]*esPq[j*esNq+m];
      }
      esE[n*esNq + m] = resE;
    }
  }
  
  //  Br = diag(nrJ.*wf);
  //  Bs = diag(nsJ.*wf);

  memory<dfloat> esBr;
  memory<dfloat> esBs;
  esBr.malloc(esNf*esNf);
  esBs.malloc(esNf*esNf);
  for(int m=0;m<esNf;++m){
    esBr[m*esNf+m] = esNrJ[m]*esWf[m];
    esBs[m*esNf+m] = esNsJ[m]*esWf[m];
  }
  
  //skew symmetric hybridized SBP operators
  //  QNr = .5*[Qr-Qr' E'*Br;
  //	    -Br*E zeros(length(wf))];
  //  QNs = .5*[Qs-Qs' E'*Bs;
  //	    -Bs*E zeros(length(wf))];
  memory<dfloat> esQNrT, esQNsT;
  esQNrT.malloc(esNp*esNp);
  esQNsT.malloc(esNp*esNp);
  
  // 0.5*(Qr-Qr')
  for(int n=0;n<esNq;++n){
    for(int m=0;m<esNq;++m){
      esQNrT[n+esNp*m] = 0.5*(esQr[n*esNq+m]-esQr[m*esNq+n]);
      esQNsT[n+esNp*m] = 0.5*(esQs[n*esNq+m]-esQs[m*esNq+n]);
    }
  }

  // E'*Br
  for(int n=0;n<esNq;++n){
    for(int m=0;m<esNf;++m){
      esQNrT[n+esNp*(m+esNq)] = 0.5*esE[m*esNq+n]*esBr[m*esNf+m];
      esQNsT[n+esNp*(m+esNq)] = 0.5*esE[m*esNq+n]*esBs[m*esNf+m];
    }
  }

  // -Br*E
  for(int n=0;n<esNf;++n){
    for(int m=0;m<esNq;++m){
      esQNrT[(n+esNq)+esNp*m] = -0.5*esE[n*esNq+m]*esBr[n*esNf+n];
      esQNsT[(n+esNq)+esNp*m] = -0.5*esE[n*esNq+m]*esBs[n*esNf+n];
    }
  }

#if 1
  memory<dfloat> esDNrT, esDNsT;
  esDNrT.malloc(esNp*esNp);
  esDNsT.malloc(esNp*esNp);
  
  for(int n=0;n<esNp;++n){
    for(int m=0;m<esNp;++m){

      esDNrT[n+m*esNp] = esQNrT[n+m*esNp] / esW[n];
      esDNsT[n+m*esNp] = esQNsT[n+m*esNp] / esW[n];

    }
  }
#endif
  
  // line up [Pq, Lf]
  memory<dfloat> esPqLfT;
  esPqLfT.malloc(mesh.Np*esNp);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      esPqLfT[n+m*mesh.Np] = esPq[n*esNq+m];
    }
    for(int m=0;m<esNf;++m){
      esPqLfT[n+(esNq+m)*mesh.Np] = esLf[n*esNf+m];
    }
  }

#if 1
  // interpolate nodes to combined quadrature
  memory<dfloat> esIqIfT;
  esIqIfT.malloc(mesh.Np*esNp);
  
  for(int n=0;n<esNp;++n){
    for(int m=0;m<mesh.Np;++m){
      if(n<esNq){
	esIqIfT[n+m*esNp] = esIq[n*mesh.Np+m];
      }else{
	esIqIfT[n+m*esNp] = esIf[(n-esNq)*mesh.Np+m];
      }
    }
  }
#endif
  
  // Fq = Vq*Pq ( vol cubature to combined quadrature through PN )
  //  dfloat *esFqT = (dfloat*) calloc(esNp*esNq, sizeof(dfloat));
  memory<dfloat> esFqT;
  esFqT.malloc(esNp*esNq);

  for(int n=0;n<esNp;++n){    
    for(int m=0;m<esNq;++m){
      dfloat resFq = 0;
      for(int j=0;j<mesh.Np;++j){
	resFq += esIqIfT[n+esNp*j]*esPq[j*esNq+m];
      }
      esFqT[n+esNp*m] = resFq;
    }
  }  
  
  memory<dfloat> esX, esY, esZ;
  esX.malloc(mesh.Nelements*esNp);
  esY.malloc(mesh.Nelements*esNp);
  esZ.malloc(mesh.Nelements*esNp);
    
  int cnt = 0;
  for(dlong e=0;e<mesh.Nelements;++e){

    dlong id = e*mesh.Nverts+0;
    
    dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh.EX[id+1];
    dfloat xe3 = mesh.EX[id+2];

    dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh.EY[id+1];
    dfloat ye3 = mesh.EY[id+2];

    for(int n=0;n<esNp;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = esR[n];
      dfloat sn = esS[n];

      /* physical coordinate of interpolation node */
      esX[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      esY[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }
  
  
  // hack in connect surface nodes here
  memory<dlong> esVmapM, esVmapP;
  esVmapM.malloc(mesh.Nelements*esNf);
  esVmapP.malloc(mesh.Nelements*esNf);
  
  for(dlong eM=0;eM<mesh.Nelements;++eM){
    for(int fM=0;fM<mesh.Nfaces;++fM){
      dlong eP = mesh.EToE[eM*mesh.Nfaces+fM];
      dlong fP = mesh.EToF[eM*mesh.Nfaces+fM];
      
      for(int n=0;n<mesh.intNfp;++n){
	int sidM  = eM*mesh.intNfp*mesh.Nfaces + fM*mesh.intNfp + n;
	int esIdM = eM*esNp + esNq + fM*mesh.intNfp + n;
	
	esVmapM[sidM] = esIdM;
	esVmapP[sidM] = esIdM; // default self connectx[

	if(!(eP==-1 || eP==eM)){
	  dfloat mindist2 = 0;
	  for(int m=0;m<mesh.intNfp;++m){
	    int sidP = eP*mesh.intNfp*mesh.Nfaces + fP*mesh.intNfp + m;

	    dfloat dist2 = pow(mesh.intx[sidP]-mesh.intx[sidM], 2) + pow(mesh.inty[sidP]-mesh.inty[sidM], 2);

	    if(dist2<mindist2 || m==0){
	      int esIdP = eP*esNp + esNq + fP*mesh.intNfp + m;
	      
	      mindist2 = dist2;
	      esVmapP[sidM] = esIdP;	      
	    }
	  }
	}
      }
    }
  }

  memory<dfloat> esMu;
  esMu.malloc(esNp*mesh.Nelements);
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<esNp;++n){
      dlong id = e*esNp+n;
      dfloat mu = meanMu;
      if(outflowMu)
	mu += 0.5*(1+tanh(10*(esX[id]-0.5*(outflowXmin+outflowXmax))/(outflowXmax-outflowXmin)))*(outflowMu-meanMu);
      esMu[id] = mu;
    }
  }
  o_esMu = platform.malloc<dfloat>(esNp*mesh.Nelements, esMu);

  o_esSurfRHS = platform.malloc<dfloat>(Nfields*mesh.Nelements*esNp);
  
  o_esX = platform.malloc<dfloat>(mesh.Nelements*esNp, esX);
  o_esY = platform.malloc<dfloat>(mesh.Nelements*esNp, esY);
  o_esZ = platform.malloc<dfloat>(mesh.Nelements*esNp, esZ);
  
  o_esVmapM = platform.malloc<dlong>(mesh.Nelements*mesh.intNfp*mesh.Nfaces, esVmapM);
  o_esVmapP = platform.malloc<dlong>(mesh.Nelements*mesh.intNfp*mesh.Nfaces, esVmapP);
  
  o_esIqT   = platform.malloc<dfloat>(mesh.Np*esNq,            esIqT);
  o_esIqfT  = platform.malloc<dfloat>(mesh.Np*esNp,            esIqfT);
  o_esQNrT  = platform.malloc<dfloat>(esNp*esNp, esDNrT); // TW: note switcheroo
  o_esQNsT  = platform.malloc<dfloat>(esNp*esNp, esDNsT);
  o_esPqT   = platform.malloc<dfloat>(esNq*mesh.Np,            esPqT);
  o_esPqLfT = platform.malloc<dfloat>(mesh.Np*esNp,     esPqLfT);
  o_esFqT   = platform.malloc<dfloat>(esNp*esNq,        esFqT);
  o_esLfT   = platform.malloc<dfloat>(mesh.Np*esNf,            esLfT);

  o_esItMT   = platform.malloc<dfloat>(esNq*mesh.Np,            esItMT);

  o_esQc      = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esQe      = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esQp      = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  
  o_esQcrecon = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);

  o_esR = platform.malloc<dfloat>(mesh.Np, mesh.r);
  o_esS = platform.malloc<dfloat>(mesh.Np, mesh.s);
  
  o_esRq = platform.malloc<dfloat>(esNq, esRq);
  o_esSq = platform.malloc<dfloat>(esNq, esSq);
  o_esWq = platform.malloc<dfloat>(esNq, esWq);
  o_esWf = platform.malloc<dfloat>(esNf, esWf);

  o_esdQedx = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esdQedy = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esdQedz = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);

  o_esDQe   = platform.malloc<dfloat>(mesh.Np*Nfields*mesh.Nelements);
  
  
  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["includes"] += DESDG "okl/esdgHelpers.h";

  kernelInfo["defines/" "p_Nverts"]= mesh.Nverts;
  
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_cubNp"]= mesh.cubNp;

  kernelInfo["defines/" "p_fluxNp"]= fluxNp;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 512/mesh.Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // entropy stable kernels
  kernelInfo["defines/" "p_esNp"] = esNp;
  kernelInfo["defines/" "p_esNq"] = esNq;
  kernelInfo["defines/" "p_esNf"] = esNf;

  // bounds guards
  int guardLevel = 0;
  settings.getSetting("GUARD LEVEL",guardLevel);
  kernelInfo["defines/" "GUARD LEVEL"]= guardLevel;
  
  // mean flow settings
  dfloat rbar, ubar, vbar, pbar;
  settings.getSetting("MEAN FLOW DENSITY", rbar);
  settings.getSetting("MEAN FLOW XVELOCITY", ubar);
  settings.getSetting("MEAN FLOW YVELOCITY", vbar);
  settings.getSetting("MEAN FLOW PRESSURE", pbar);
  kernelInfo["defines/" "p_rbar"] = rbar;
  kernelInfo["defines/" "p_ubar"] = ubar;
  kernelInfo["defines/" "p_vbar"] = vbar;
  kernelInfo["defines/" "p_pbar"] = pbar;
  
  //  std::cout << kernelInfo << std::endl;
  
  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

#if 0
  // kernels from volume file
  sprintf(fileName, DESDG "/okl/esdgCubatureVolume%s.okl", suffix);
  sprintf(kernelName, "esdgCubatureVolume%s", suffix);
  
  cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
						     kernelInfo);
  // kernels from surface file
  sprintf(fileName, DESDG "/okl/esdgCubatureSurface%s.okl", suffix);
  sprintf(kernelName, "esdgCubatureSurface%s", suffix);
  
  cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
					     kernelInfo);

  if(mesh.dim==2){
    sprintf(fileName, DESDG "/okl/esdgCubatureInitialCondition2D.okl");
    sprintf(kernelName, "esdgCubatureInitialCondition2D");
  }
  
  if(mesh.dim==3){
    sprintf(fileName, DESDG "/okl/esdgCubatureInitialCondition3D.okl");
    sprintf(kernelName, "esdgCubatureInitialCondition3D");
  }

  cubatureInitialConditionKernel = platform.buildKernel(fileName, kernelName,
						      kernelInfo);


  // cubature surface interpolation kernel
  sprintf(fileName, DESDG "/okl/esdgCubatureSurfaceInterpolation%s.okl", suffix);
  sprintf(kernelName, "esdgCubatureSurfaceInterpolation%s", suffix);
  
  cubatureSurfaceInterpolationKernel = platform.buildKernel(fileName, kernelName,
							  kernelInfo);

#endif
  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/libs/mesh/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  MassMatrixKernel = platform.buildKernel(fileName, kernelName,
					kernelInfo);


  sprintf(fileName, DESDG "/okl/esdgError%s.okl", suffix);
  sprintf(kernelName, "esdgError%s", suffix);
  
  errorKernel =  platform.buildKernel(fileName, kernelName,
				    kernelInfo);

  
  // vorticity calculation
  sprintf(fileName, DESDG "/okl/esdgVorticity%s.okl", suffix);
  sprintf(kernelName, "esdgVorticity%s", suffix);

  vorticityKernel = platform.buildKernel(fileName, kernelName,
				       kernelInfo);

  if (mesh.dim==2) {

    sprintf(fileName, DESDG "/okl/esdgDGVorticity%s.okl", suffix);
    sprintf(kernelName, "esdgDGVorticity%s", suffix);
    
    dgVorticityKernel = platform.buildKernel(fileName, kernelName,
					   kernelInfo);

    printf("NODAL 2D INITIAL CONDITION\n");    
    sprintf(fileName, DESDG "/okl/esdgInitialCondition2D.okl");
    sprintf(kernelName, "esdgInitialCondition2D");
  } else {
    printf("NODAL 3D INITIAL CONDITION\n");
    sprintf(fileName, DESDG "/okl/esdgInitialCondition3D.okl");
    sprintf(kernelName, "esdgInitialCondition3D");
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
					      kernelInfo);

  
  sprintf(fileName, DESDG "/okl/esdg%s.okl", suffix);

  sprintf(kernelName, "esdgInterpolate%s", suffix);

  esInterpolateKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  sprintf(kernelName, "esdgIntegrateEntropyChange%s", suffix);

  esIntegrateEntropyChangeKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  sprintf(kernelName, "esdgIntegrateEntropy%s", suffix);

  esIntegrateEntropyKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  
  sprintf(kernelName, "esdgVolume%s", suffix);
  
  esVolumeKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  sprintf(kernelName, "esdgSurface%s", suffix);
  
  esSurfaceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  
  sprintf(kernelName, "esdgVolumeGradient%s", suffix);
  
  esVolumeGradientKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  sprintf(kernelName, "esdgSurfaceGradient%s", suffix);
  
  esSurfaceGradientKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  sprintf(kernelName, "esdgVolumeDivergence%s", suffix);
  
  esVolumeDivergenceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  sprintf(kernelName, "esdgSurfaceDivergence%s", suffix);
  
  esSurfaceDivergenceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  sprintf(kernelName, "esdgDiffusionFluxes%s", suffix);
  
  esDiffusionFluxesKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);
  
}

