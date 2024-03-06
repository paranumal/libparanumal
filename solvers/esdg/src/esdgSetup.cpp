/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

void esdg_t::Setup(platform_t& _platform, mesh_t& _mesh,
                  esdgSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //get physical paramters
  settings.getSetting("VISCOSITY", mu);
  settings.getSetting("GAMMA", gamma);

  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  isothermal = (settings.compareSetting("ISOTHERMAL", "TRUE")) ? 1:0;

  //setup cubature
  mesh.CubatureSetup();
  mesh.CubaturePhysicalNodes();

  Nfields   = mesh.dim + 2;
  Ngrads = mesh.dim*mesh.dim;

  // build element centers
  memory<dfloat> cx(mesh.Nelements, 0);
  memory<dfloat> cy(mesh.Nelements, 0);
  memory<dfloat> cz(mesh.Nelements, 0);

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

  o_cx = platform.malloc<dfloat>(cx);
  o_cy = platform.malloc<dfloat>(cy);
  o_cz = platform.malloc<dfloat>(cz);
  
  // Build ESDG specific connectivities
  memory<int> NVToE(mesh.Nnodes, 0);
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int v=0;v<mesh.Nverts;++v){
      ++NVToE[mesh.EToV[e*mesh.Nverts+v]];
    }
  }
  int maxNVToE = 0;
  for(dlong n=0;n<mesh.Nnodes;++n){
    maxNVToE = std::max(maxNVToE, NVToE[n]);
    NVToE[n] = 0;
  }

  memory<dlong> VToE(mesh.Nnodes*maxNVToE);
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int v=0;v<mesh.Nverts;++v){
      dlong id = mesh.EToV[e*mesh.Nverts+v];
      VToE[id*maxNVToE+NVToE[id]] = e;
      ++NVToE[id];
    }
  }

  // list of local vertices to elements (via vertex connectivity)
  memory<dlong> LVToE(mesh.Nelements*mesh.Nverts*maxNVToE, 0);
  memory<dlong> NLVToE(mesh.Nelements*mesh.Nverts, 0);
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int v=0;v<mesh.Nverts;++v){
      dlong id = mesh.EToV[e*mesh.Nverts+v];
      for(int n=0;n<NVToE[id];++n){
	LVToE[e*mesh.Nverts*maxNVToE + v*maxNVToE + n] = VToE[id*maxNVToE+n];
      }
      NLVToE[e*mesh.Nverts+v] = NVToE[id];
    }
  }
  
#if 0
  #include "../../libs/mesh/include/meshDefines2D.h"
  dfloat minh = 1e9, maxh = 0;
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat Je = mesh.vgeo[e*mesh.Nvgeo+JID];
    minh = mymin(minh, sqrt(Je));
    maxh = mymax(maxh, sqrt(Je));
  }
  printf("h in range [%g, %g]\n", minh, maxh);
#endif
  
  
  // sort layout later
  maxLVToE = maxNVToE;
  o_LVToE = platform.malloc<dlong>(LVToE);
  o_NLVToE = platform.malloc<dlong>(NLVToE);
  
  memory<int> NEToVToE(mesh.Nelements, 0);
  maxEToVToE = maxNVToE*mesh.Nverts;

  memory<dlong> EToVToE(maxEToVToE*mesh.Nelements, 0);
  
  // use hlongs to get away with mysort
  memory<hlong> tmpEToVToE(maxEToVToE, 0);

  for(dlong e=0;e<mesh.Nelements;++e){

    // first count
    int Nconns = 0;
    for(int v=0;v<mesh.Nverts;++v){
      dlong id = mesh.EToV[e*mesh.Nverts+v];
      for(int c=0;c<NVToE[id];++c){
	tmpEToVToE[Nconns++] = VToE[id*maxNVToE+c];
      }
    }

    // sort temporary EToE for this element
    //    mysort(tmpEToVToE, Nconns, "ascend");
    std::sort(tmpEToVToE.ptr(), tmpEToVToE.ptr()+Nconns, std::less<hlong>());

    // remove dups
    int finalNconns = 1;
    for(int n=1;n<Nconns;++n){
      if(tmpEToVToE[n]!=tmpEToVToE[finalNconns-1]){
	tmpEToVToE[finalNconns++] = tmpEToVToE[n];
      }
    }

    // record number of EToVToE connections
    NEToVToE[e] = finalNconns;
    for(int n=0;n<finalNconns;++n){
      EToVToE[e*maxEToVToE+n] = tmpEToVToE[n];
    }
  }

  o_NEToVToE = platform.malloc(NEToVToE);
  o_EToVToE  = platform.malloc(EToVToE);
  
  // port EToE to device
  o_EToE = platform.malloc(mesh.EToE);
  o_EToB = platform.malloc(mesh.EToB);


  // face cubature
  // build interpolation matrix: nodes to surface nodes
  // HACKITY-HACK (not saved anywhere)
  memory<dfloat> ir(mesh.intNfp*mesh.Nfaces, 0);
  memory<dfloat> is(mesh.intNfp*mesh.Nfaces, 0);
  memory<dfloat> it(mesh.intNfp*mesh.Nfaces, 0);
  memory<dfloat> iw(mesh.intNfp*mesh.Nfaces, 0);
  
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

  }

  if(mesh.dim==3){
    // build interpolation matrix: nodes to surface nodes
    for(int n=0;n<mesh.intNfp;++n){
      ir[0*mesh.intNfp + n] =  mesh.intr[n];
      ir[1*mesh.intNfp + n] =  mesh.intr[n];
      ir[2*mesh.intNfp + n] =  mesh.intr[n];
      ir[3*mesh.intNfp + n] = -1.0;
      
      is[0*mesh.intNfp + n] =  mesh.ints[n];
      is[1*mesh.intNfp + n] = -1.0;
      is[2*mesh.intNfp + n] =  mesh.ints[n];
      is[3*mesh.intNfp + n] =  mesh.intr[n];
      
      it[0*mesh.intNfp + n] = -1.0;
      it[1*mesh.intNfp + n] =  mesh.ints[n];
      it[2*mesh.intNfp + n] = -(1.0 + mesh.intr[n] + mesh.ints[n]);
      it[3*mesh.intNfp + n] =  mesh.ints[n];
      
      iw[0*mesh.intNfp + n] =  mesh.intw[n];
      iw[1*mesh.intNfp + n] =  mesh.intw[n];
      iw[2*mesh.intNfp + n] =  mesh.intw[n];
      iw[3*mesh.intNfp + n] =  mesh.intw[n];
    }

  }

  
  // entropy stable infrastructure following On discretely entropy conservative and entropy stable discontinuous Galerkin methods, Jesse Chan
  
  int esNq = mesh.cubNp;

  memory<dfloat> esRq(esNq, 0);
  memory<dfloat> esSq(esNq, 0);
  memory<dfloat> esWq(esNq, 0);

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
  memory<dfloat> esRf(esNf, 0);
  memory<dfloat> esSf(esNf, 0);
  memory<dfloat> esWf(esNf, 0);

  for(int n=0;n<esNf;++n){
    esRf[n] = ir[n];
    esSf[n] = is[n];
    esWf[n] = iw[n];
  }

  memory<dfloat> esR(esNp, 0);
  memory<dfloat> esS(esNp, 0);
  memory<dfloat> esW(esNp, 0);
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
  memory<dfloat> esIq(mesh.Np*esNp, 0);
  memory<dfloat> esIqT(mesh.Np*esNp, 0);
  memory<dfloat> esIf(mesh.Np*esNf, 0);

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
  memory<dfloat> esIqfT(mesh.Np*esNp, 0);
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
  memory<dfloat> esMM(mesh.Np*esNq, 0);
  memory<dfloat> esInvMM(mesh.Np*esNq, 0);
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

  memory<dfloat> esPq (mesh.Np*esNq, 0);
  memory<dfloat> esPqT(mesh.Np*esNq, 0);
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

  memory<dfloat> esItMT (mesh.Np*esNq, 0);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      dfloat resItM = 0;

      esItMT[m+esNq*n] = esWq[m]*esIq[m*mesh.Np+n];
    }
  }


  
  // Lf = M\(Vf'*diag(wf));
  memory<dfloat> esLf(mesh.Np*esNf, 0);
  memory<dfloat> esLfT(mesh.Np*esNf, 0);
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

  memory<dfloat> esIqfLfT(esNp*esNf, 0);
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

  o_esIqfLfT = platform.malloc<dfloat>(esIqfLfT);


  memory<dfloat> esIqfDrPqT(esNp*esNq, 0);
  memory<dfloat> esIqfDsPqT(esNp*esNq, 0);
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

  o_esIqfDrPqT = platform.malloc<dfloat>(esIqfDrPqT);
  o_esIqfDsPqT = platform.malloc<dfloat>(esIqfDsPqT);


  // scaled reference normals
  // z = zeros(size(rq1D));
  // e = ones(size(rq1D));
  // nrJ = [-z; e; -e];
  // nsJ = [-e; e; -z];

  memory<dfloat> esNrJ(esNf, 0);
  memory<dfloat> esNsJ(esNf, 0);
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
  memory<dfloat> esQrModal(mesh.Np*mesh.Np, 0);
  memory<dfloat> esQsModal(mesh.Np*mesh.Np, 0);

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

  memory<dfloat> esQr(esNq*esNq, 0);
  memory<dfloat> esQs(esNq*esNq, 0);

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
  memory<dfloat> esE(esNf*esNq, 0);
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

  memory<dfloat> esBr(esNf*esNf, 0);
  memory<dfloat> esBs(esNf*esNf, 0);
  for(int m=0;m<esNf;++m){
    esBr[m*esNf+m] = esNrJ[m]*esWf[m];
    esBs[m*esNf+m] = esNsJ[m]*esWf[m];
  }
  
  //skew symmetric hybridized SBP operators
  //  QNr = .5*[Qr-Qr' E'*Br;
  //	    -Br*E zeros(length(wf))];
  //  QNs = .5*[Qs-Qs' E'*Bs;
  //	    -Bs*E zeros(length(wf))];
  memory<dfloat> esQNrT(esNp*esNp, 0);
  memory<dfloat> esQNsT(esNp*esNp, 0);

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

  memory<dfloat> esDNrT(esNp*esNp, 0);
  memory<dfloat> esDNsT(esNp*esNp, 0);
  
  for(int n=0;n<esNp;++n){
    for(int m=0;m<esNp;++m){

      esDNrT[n+m*esNp] = esQNrT[n+m*esNp] / esW[n];
      esDNsT[n+m*esNp] = esQNsT[n+m*esNp] / esW[n];

    }
  }
  
  // line up [Pq, Lf]
  memory<dfloat> esPqLfT(mesh.Np*esNp, 0);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      esPqLfT[n+m*mesh.Np] = esPq[n*esNq+m];
    }
    for(int m=0;m<esNf;++m){
      esPqLfT[n+(esNq+m)*mesh.Np] = esLf[n*esNf+m];
    }
  }

  // interpolate nodes to combined quadrature
  memory<dfloat> esIqIfT(esNp*mesh.Np, 0);
  for(int n=0;n<esNp;++n){
    for(int m=0;m<mesh.Np;++m){
      if(n<esNq){
	esIqIfT[n+m*esNp] = esIq[n*mesh.Np+m];
      }else{
	esIqIfT[n+m*esNp] = esIf[(n-esNq)*mesh.Np+m];
      }
    }
  }
  
  // Fq = Vq*Pq ( vol cubature to combined quadrature through PN )
  memory<dfloat> esFqT(esNp*esNq, 0);
  for(int n=0;n<esNp;++n){    
    for(int m=0;m<esNq;++m){
      dfloat resFq = 0;
      for(int j=0;j<mesh.Np;++j){
	resFq += esIqIfT[n+esNp*j]*esPq[j*esNq+m];
      }
      esFqT[n+esNp*m] = resFq;
    }
  }  

  /* build artificial relaxation operator */
  memory<dfloat> esVq (esNq*mesh.Np, 0);
  memory<dfloat> esVWB(mesh.Np*mesh.Np, 0);
  
  mesh.VandermondeTri2D(mesh.N, esRq,   esSq, esVq);
  mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, esVWB);

  memory<dfloat> esRelaxD  (mesh.Np, 0);

  int sk = 0;
  // borrowed directly from VDMTri2D
  for(int i=0; i<mesh.N+1; i++){
    for(int j=0; j<mesh.N+1-i; j++){
      // set relaxation rate for each mode here
      int cut = 0;
      if(i+j>cut){
	esRelaxD[sk] = -2*(i+j); // pow((dfloat)(i+j-cut)/(dfloat)(mesh.N-cut), 4);
      }else{
	esRelaxD[sk] = 0;
      }
      printf("esRelaxD[%d (%d,%d)] = %g\n", sk, i, j, esRelaxD[sk]);
      ++sk;
    }
  }

  memory<dfloat> esRelaxOpT(esNq*mesh.Np, 0);
  for(int n=0;n<mesh.Np;++n){
    for(int m=0;m<esNq;++m){
      dfloat res = 0;
      for(int i=0;i<mesh.Np;++i){
	res += esVWB[n*mesh.Np+i]*esRelaxD[i]*esVq[m*mesh.Np+i]*esWq[m]; // (VWB*D*Vq'*wq)
      }
      esRelaxOpT[n + m*mesh.Np] = res;
    }
  }

  o_esRelaxOpT = platform.malloc<dfloat>(esRelaxOpT);
  
  memory<dfloat> esX(mesh.Nelements*esNp, 0);
  memory<dfloat> esY(mesh.Nelements*esNp, 0);
  memory<dfloat> esZ(mesh.Nelements*esNp, 0);

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
  memory<dlong> esVmapM(mesh.Nelements*esNf,0);
  memory<dlong> esVmapP(mesh.Nelements*esNf,0);
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

  dfloat meanMu;
  settings.getSetting("LAME MU", meanMu);
  
  // skip
  memory<dfloat> esMu(esNp*mesh.Nelements,0);
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<esNp;++n){
      dlong id = e*esNp+n;
      dfloat mu = meanMu;
      //      if(outflowMu)
      //	mu += 0.5*(1+tanh(10*(esX[id]-0.5*(outflowXmin+outflowXmax))/(outflowXmax-outflowXmin)))*(outflowMu-meanMu);
      esMu[id] = mu;
    }
  }
  o_esMu = platform.malloc<dfloat>(esMu);

  o_esX = platform.malloc<dfloat>(esX);
  o_esY = platform.malloc<dfloat>(esY);
  o_esZ = platform.malloc<dfloat>(esZ);

  o_esVmapM = platform.malloc<dlong>(esVmapM);
  o_esVmapP = platform.malloc<dlong>(esVmapP);
  	    
  o_esIqT   = platform.malloc<dfloat>(esIqT);
  o_esIqfT  = platform.malloc<dfloat>(esIqfT);
  o_esQNrT  = platform.malloc<dfloat>(esDNrT); // TW: note switcheroo
  o_esQNsT  = platform.malloc<dfloat>(esDNsT);
  o_esPqT   = platform.malloc<dfloat>(esPqT);
  o_esPqLfT = platform.malloc<dfloat>(esPqLfT);
  o_esFqT   = platform.malloc<dfloat>(esFqT);
  o_esLfT   = platform.malloc<dfloat>(esLfT);
	    
  o_esItMT  = platform.malloc<dfloat>(esItMT);
  	    
  o_esR     = platform.malloc<dfloat>(mesh.r);
  o_esS     = platform.malloc<dfloat>(mesh.s);
  	    
  o_esRq    = platform.malloc<dfloat>(esRq);
  o_esSq    = platform.malloc<dfloat>(esSq);
  o_esWq    = platform.malloc<dfloat>(esWq);
  o_esWf    = platform.malloc<dfloat>(esWf);

  o_esQc    = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esQe    = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esQp    = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  
  o_esQcrecon = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  
  o_esdQedx = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esdQedy = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
  o_esdQedz = platform.malloc<dfloat>(esNp*Nfields*mesh.Nelements);
	    
  o_esDQe   = platform.malloc<dfloat>(mesh.Np*Nfields*mesh.Nelements);

  o_esSurfRHS = platform.malloc<dfloat>(Nfields*mesh.Nelements*esNp);


  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "max"});

  /*setup trace halo exchange */
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  dlong NlocalFields = mesh.Nelements*mesh.Np*Nfields;
  dlong NhaloFields  = mesh.totalHaloPairs*mesh.Np*Nfields;

  // compute samples of q at interpolation nodes
  q.malloc(NlocalFields+NhaloFields);
  o_q = platform.malloc<dfloat>(NlocalFields+NhaloFields);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["includes"] += DESDG "okl/esdgHelpers.h";

  std::cout << kernelInfo << std::endl;
  
  kernelInfo["defines/" "p_Nverts"]= mesh.Nverts;
  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  kernelInfo["defines/" "p_cubNp"]= mesh.cubNp;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  if (cubature) {
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1, blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1, blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }

  // entropy stable kernels
  kernelInfo["defines/" "p_esNp"] = esNp;
  kernelInfo["defines/" "p_esNq"] = esNq;
  kernelInfo["defines/" "p_esNf"] = esNf;

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
  
  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DESDG "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName = oklFilePrefix + "esdg" + suffix + oklFileSuffix;

  kernelName = "esdgInterpolate" + suffix;

  esInterpolateKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "esdgIntegrateEntropyChange" + suffix;

  esIntegrateEntropyChangeKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "esdgIntegrateEntropy" + suffix;

  esIntegrateEntropyKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  
  kernelName = "esdgVolume" +  suffix;
  
  esVolumeKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "esdgSurface" + suffix;
  
  esSurfaceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  
  kernelName = "esdgVolumeGradient" + suffix;
  
  esVolumeGradientKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "esdgSurfaceGradient" + suffix;
  
  esSurfaceGradientKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "esdgVolumeDivergence" + suffix;
  
  esVolumeDivergenceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "esdgSurfaceDivergence" + suffix;
  
  esSurfaceDivergenceKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName =  "esdgDiffusionFluxes" + suffix;
  
  esDiffusionFluxesKernel =
    platform.buildKernel(fileName, kernelName, kernelInfo);

  // vorticity calculation
  fileName   = oklFilePrefix + "esdgVorticity" + suffix + oklFileSuffix;
  kernelName = "esdgVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                     kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "esdgInitialCondition2D" + oklFileSuffix;
    kernelName = "esdgInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "esdgInitialCondition3D" + oklFileSuffix;
    kernelName = "esdgInitialCondition3D";
  }
  
  initialConditionKernel = platform.buildKernel(fileName, kernelName,
						kernelInfo);

#if 0
  fileName   = oklFilePrefix + "esdgMaxWaveSpeed" + suffix + oklFileSuffix;
  if (isothermal) {
    kernelName = "esdgIsothermalMaxWaveSpeed" + suffix;
  } else {
    kernelName = "esdgMaxWaveSpeed" + suffix;
  }

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName,
                                            kernelInfo);
#endif
}
