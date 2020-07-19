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

#include "subcell.hpp"
#include "subcell2D.hpp"
// #include "mesh3D.hpp"

// subcellTri2D::subcellTri2D(mesh_t& _mesh, settings_t& _settings):
//    subcell_t(_mesh, _settings) {
//    _settings.getSetting("SUBCELL NUMBER", N);
//     // Nverts = mesh.Nverts; // Subcell could be in different topology but not now
//     // Nfaces = mesh.Nfaces; 
//     // // Currently equispaced triangle
//     // Nsubcells = N*N;  // 
//     // // Number of nodes
//     // Np = 0.5*(N+1)*(N+2);
//     //
//    }

subcellTri2D::subcellTri2D(solver_t& _solver):
   subcell_t(_solver) {
   settings.getSetting("SUBCELL NUMBER", N);
    // Nverts = mesh.Nverts; // Subcell could be in different topology but not now
    // Nfaces = mesh.Nfaces; 
    // // Currently equispaced triangle
    // Nsubcells = N*N;  // 
    // // Number of nodes
    // Np = 0.5*(N+1)*(N+2);
    //
   }

void subcellTri2D::SetupDetector(){

  // Set mode map info for skyline 
  Nmodes   = mesh.Np; 
  Nmodes1D = mesh.N+1; 

  ModeMap  = (int *) malloc(Nmodes*sizeof(int));
  ModeInfoTri2D(mesh.N, ModeMap); 
  //

  LSF = (dfloat *) malloc(mesh.N*sizeof(dfloat));
  LeastSquaresFit(mesh.N, LSF);

  BLD = (dfloat *) malloc((mesh.N+1)*sizeof(dfloat));
  BaseLineDecay(mesh.N, BLD);

  // Get inverse Vandermonde Matrix 
  dfloat *invV = (dfloat *) malloc(mesh.Np*mesh.Np*sizeof(dfloat));
  mesh.VandermondeTri2D(mesh.N, mesh.Np, mesh.r, mesh.s, invV);
  matrixInverse(mesh.Np, invV); 
  // Create transpose
  invVT = (dfloat *) malloc(mesh.Np*mesh.Np*sizeof(dfloat));
  for(int m=0; m<mesh.Np; m++)
    for(int n=0; n<mesh.Np; n++)
      invVT[n+m*mesh.Np]  = invV[m+n*mesh.Np]; 
  
  // Initialize Element List on host
  ElementList = (dlong *) calloc(mesh.Nelements, sizeof(dlong)); 

}


void subcellTri2D::OccaSetup(){
 // occa::properties kernelInfo = props; //copy base properties
 occa::properties &kernelInfo = props; 
 // Detector Related
 o_ElementList = device.malloc(mesh.Nelements*sizeof(dlong), ElementList);
 // o_ElementList = device.malloc(mesh.Nelements*sizeof(dfloat), ElementList);
 o_invVT       = device.malloc(mesh.Np*mesh.Np*sizeof(dfloat), invVT); 
 o_ModMap      = device.malloc(mesh.Np*sizeof(int), ModeMap); 
 o_LSF         = device.malloc(mesh.N*sizeof(dfloat), LSF); 
 o_BLD         = device.malloc((mesh.N+1)*sizeof(dfloat), BLD); 

 // Projection Matrix
 dfloat *PMT = (dfloat *) malloc(mesh.Np*Nsubcells*sizeof(dfloat)); 
 for(int n=0; n<Nsubcells; n++){
  for(int m=0; m<mesh.Np; m++){
    PMT[m*Nsubcells +n] = PM[n*mesh.Np + m]; 
  }
 }
 o_PMT  = device.malloc(mesh.Np*Nsubcells*sizeof(dfloat), PMT);
 free(PMT); 

// Recontruction Matrix
 dfloat *RMT = (dfloat *) malloc(mesh.Np*Nsubcells*sizeof(dfloat)); 
 for(int n=0; n<mesh.Np; n++){
  for(int m=0; m<Nsubcells; m++){
    RMT[m*mesh.Np +n] = RM[n*Nsubcells + m]; 
  }
 }

 o_RMT  = device.malloc(mesh.Np*Nsubcells*sizeof(dfloat), RMT);
 free(RMT);

 // Face projection matrix
 dfloat *PFMT = (dfloat*) malloc(N*mesh.Nfp*sizeof(dfloat));
 for(int n=0; n<N; n++){
  for(int m=0; m<mesh.Nfp; m++){
    PFMT[m*N + n] = PFM[n*mesh.Nfp + m]; 
  }
 }
o_PFMT  = device.malloc(N*mesh.Nfp*sizeof(dfloat), PFMT);
free(PFMT);

 
 // Face reconstruct matrix 
 dfloat *RFMT = (dfloat*) malloc(N*mesh.Nfp*sizeof(dfloat));
 for(int n=0; n<mesh.Nfp; n++){
  for(int m=0; m<N; m++){
    RFMT[m*mesh.Nfp+n] = RFM[n*N + m]; 
  }
 }
o_RFMT  = device.malloc(N*mesh.Nfp*sizeof(dfloat), RFMT);
free(RFMT);


// Face reconstruct matrix 
 dfloat *SLIFTT = (dfloat*) malloc(mesh.Np*N*mesh.Nfaces*sizeof(dfloat));
 for(int n=0; n<mesh.Np; n++){
  for(int m=0; m<(N*Nfaces); m++){
    SLIFTT[m*mesh.Np+n] = SLIFT[n*(N*Nfaces) + m]; 
  }
 }
o_SLIFTT  = device.malloc(mesh.Np*N*Nfaces*sizeof(dfloat), SLIFTT);
free(SLIFTT);


 
 const dlong Nlocal = mesh.Nelements*Nsubcells; 
 const dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*Nsubcells;  
 o_vgeo = device.malloc(Ntotal*Nvgeo*sizeof(dfloat), vgeo); 
 o_sgeo = device.malloc(Ntotal*Nfaces*Nsgeo*sizeof(dfloat), sgeo); 


 o_EToE = device.malloc(mesh.Nelements*mesh.Nfaces*sizeof(dlong), mesh.EToE);
 // Connectivity
 o_ielist = device.malloc(Nint*sizeof(int), ielist);
 o_eelist = device.malloc(Next*sizeof(int), eelist);

 // o_iemapP = device.malloc(Nint*Nfaces*sizeof(int), iemapP);
 // o_ifmapP = device.malloc(Nint*Nfaces*sizeof(int), ifmapP);

 // o_eemapP = device.malloc(Next*Nfaces*mesh.Nelements*sizeof(dlong), eemapP); 
 // o_efmapP = device.malloc(Next*Nfaces*mesh.Nelements*sizeof(dlong), efmapP); 
 o_emapP = device.malloc(Nsubcells*Nfaces*mesh.Nelements*sizeof(dlong), emapP); 
 o_fmapP = device.malloc(Nsubcells*Nfaces*mesh.Nelements*sizeof(dlong), fmapP); 

 o_mFToE = device.malloc(Nfaces*N*sizeof(int), mFToE);
 o_mFToF = device.malloc(Nfaces*N*sizeof(int), mFToF);
 o_mDGID = device.malloc(Nfaces*N*sizeof(int), mDGID);

  // Some Defines
  kernelInfo["defines/" "s_Nvgeo"] = Nvgeo;  
  kernelInfo["defines/" "s_Nsgeo"] = Nsgeo;  
  kernelInfo["defines/" "s_Ncells"]= Nsubcells;  
  kernelInfo["defines/" "s_N"]     = N; 
  kernelInfo["defines/" "s_Nint"]  = Nint; 
  kernelInfo["defines/" "s_Next"]  = Next; 
  

  kernelInfo["defines/" "s_DGDG_TYPE"]  = (int)DGDG_TYPE; 
  kernelInfo["defines/" "s_DGFV_TYPE"]  = (int)DGFV_TYPE; 
  kernelInfo["defines/" "s_FVFV_TYPE"]  = (int)FVFV_TYPE; 

  
  int maxNodes = std::max(Nsubcells, mesh.Np); 
  kernelInfo["defines/" "s_maxNodes"] = maxNodes;  

  kernelInfo["defines/" "s_NfacesNfp"] = Nfaces*N;   

  int Nblock = 1; //1024/ maxNodes; 
  kernelInfo["defines/" "s_Nblocks"] = Nblock; 

  int maxNodesS = std::max(Nfaces*mesh.Nfp, Nfaces*N); 
  kernelInfo["defines/" "s_maxNodesS"] = maxNodesS;  
  
  int NblockS = 1;  // 512/ maxNodesS; 
  kernelInfo["defines/" "s_NblocksS"] = NblockS; 

  //
  kernelInfo["defines/" "s_CXID"]= CXID;
  kernelInfo["defines/" "s_CYID"]= CYID;
  kernelInfo["defines/" "s_IVID"]= IVID;

  kernelInfo["defines/" "s_FXID"]= FXID;
  kernelInfo["defines/" "s_FYID"]= FYID;
  kernelInfo["defines/" "s_NXID"]= NXID;
  kernelInfo["defines/" "s_NYID"]= NYID;
  kernelInfo["defines/" "s_SAID"]= SAID;
  // kernelInfo["defines/" "s_WFID"]= WFID;
  kernelInfo["defines/" "s_BCID"]= BCID;

  
  // if(mesh.rank==0)
  //   printf("Building Kernels \n");
  // if(Nfields ==1)
  // skylineKernel = buildKernel(device, SUBCELL_DIR "/okl/"
  //                                       "detectorTri2D.okl",
  //                                       "skylineTri2D",
  //                                       kernelInfo, comm);

}


void subcellTri2D::CreateMinorGrid(){

if(settings.compareSetting("SUBCELL MINOR GRID","EQUISPACED")){
  // Using triangle does not have to be!!!!
  Nverts = mesh.Nverts; // Subcell could be in different topology but not now
  Nfaces = mesh.Nfaces; 
  // Currently equispaced triangle
  Nsubcells = N*N;  // 
  // Number of nodes in this triangulation
  Np = 0.5*(N+1)*(N+2);
  // 
  NfaceVertices = mesh.NfaceVertices; 
  faceVertices  = mesh.faceVertices; 

  vr = (dfloat *)malloc(Np*sizeof(dfloat));
  vs = (dfloat *)malloc(Np*sizeof(dfloat));
  // Currently very simple tesselation
  mesh.EquispacedNodesTri2D(N, vr, vs); 

  // Create EToV
  mEToV = (int*) malloc(Nsubcells*Nverts*sizeof(int));
  // mFToE = (int*) malloc(Nfaces*N*sizeof(int)); 
  // mFToF = (int*) malloc(Nfaces*N*sizeof(int)); 

  EquispacedEToVTri2D(N, mEToV);

  // Create Local Face to Element and Face To Face; 
}else{
  printf("This minor grid method has not been implemented yet \n ");
}


// Create basic local data i.e. center, face centers
cr = (dfloat *) malloc(Nsubcells*sizeof(dfloat));
cs = (dfloat *) malloc(Nsubcells*sizeof(dfloat));
// 
fr = (dfloat *) malloc(Nsubcells*Nfaces*sizeof(dfloat));
fs = (dfloat *) malloc(Nsubcells*Nfaces*sizeof(dfloat));

mJ = (dfloat *)malloc(Nsubcells*sizeof(dfloat));

for(int s=0; s<Nsubcells; s++){
  dfloat tmpx = 0.0, tmpy = 0.0;
  // center 
  for(int v=0; v<Nverts; v++){
    const int vid = mEToV[s*Nverts+v];
    tmpx += vr[vid]; 
    tmpy += vs[vid];
  }

  //
  cr[s] = tmpx/Nverts;  
  cs[s] = tmpy/Nverts;

  for(int f=0; f<Nfaces; f++){
    tmpx = 0.0, tmpy = 0.0; 
    for(int n=0;n<NfaceVertices;++n){
      const int vid = s*Nverts + faceVertices[f*NfaceVertices+n];
      tmpx += vr[mEToV[vid]]/NfaceVertices; 
      tmpy += vs[mEToV[vid]]/NfaceVertices; 
    }
    //
    fr[s*Nfaces+f] = tmpx;
    fs[s*Nfaces+f] = tmpy;
  } 

  const int v1 = mEToV[s*Nverts+0]; 
  const int v2 = mEToV[s*Nverts+1]; 
  const int v3 = mEToV[s*Nverts+2]; 

  const dfloat xv1 = vr[v1];    
  const dfloat yv1 = vs[v1];

  const dfloat xv2 = vr[v2];    
  const dfloat yv2 = vs[v2];    
  
  const dfloat xv3 = vr[v3];    
  const dfloat yv3 = vs[v3];    
  // A_subcell/A_reference
  mJ[s] = 0.5*((xv2-xv1)*(yv3-yv1) - (xv3-xv1)*(yv2-yv1))/2.0;
  // printf("J = %.12e", mJ[s]);
}

// Required for mixed element lifting
mDGID = (int *)malloc(N*Nfaces*sizeof(int)); 
for(int n=0; n<Nfaces*N; n++){
 const int face        = n/N;  
 const int lid         = n%N;  
 if(lid<mesh.Nfp){
  const int dgid = face*N + lid; 
  mDGID[dgid]  = face*mesh.Nfp + lid; 
 }
}

#if 0
  for(int e=0; e<Nsubcells; e++)
    printf("%d %d %d \n", EToV[e*Nverts + 0]+1,EToV[e*Nverts + 1]+1,EToV[e*Nverts + 2]+1 );

  for(int f=0; f<mesh.Nfaces; f++)
    for(int n=0; n<mesh.Nfp; n++)
      printf("f: %d n= %d \n", f, mesh.faceNodes[f*mesh.Nfp + n]);

  rc = (dfloat *)malloc(Nsubcells*sizeof(dfloat));
  sc = (dfloat *)malloc(Nsubcells*sizeof(dfloat));

  for(int e=0; e<Nsubcells;e++){
    dfloat vr = 0.0, vs = 0.0; 
    for(int n=0; n<Nverts; n++){
     vr += sr[mEToV[e*Nverts + n]]; 
     vs += ss[mEToV[e*Nverts + n]]; 
    }
    rc[e] = vr/Nverts; sc[e] = vs/Nverts; 
  }
#endif
}


void subcellTri2D::GeometricFactors(){
  Nvgeo = 3; 
  Nsgeo = 6; 
  vgeo = (dfloat*) calloc((mesh.Nelements+mesh.totalHaloPairs)*Nsubcells*Nvgeo, 
                           sizeof(dfloat));

  // Compute Volume Geometric Facors
  for(dlong e=0; e<mesh.Nelements; e++){
    dlong id = e*Nverts;
    const dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
    const dfloat xe2 = mesh.EX[id+1];
    const dfloat xe3 = mesh.EX[id+2];

    const dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
    const dfloat ye2 = mesh.EY[id+1];
    const dfloat ye3 = mesh.EY[id+2];

    for(int s= 0; s<Nsubcells; s++){
      const dlong elm = e*Nsubcells + s;  
      // Compute Subcell Centers
      const dfloat rn = cr[s]; 
      const dfloat sn = cs[s]; 
      // Cell centers
      vgeo[Nvgeo*elm + CXID] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      vgeo[Nvgeo*elm + CYID] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      //
       // Get all vertices
      const int sv1 = mEToV[s*Nverts + 0]; 
      const int sv2 = mEToV[s*Nverts + 1]; 
      const int sv3 = mEToV[s*Nverts + 2];

     const dfloat rn1 = vr[sv1];    
     const dfloat sn1 = vs[sv1];    
     const dfloat rn2 = vr[sv2];    
     const dfloat sn2 = vs[sv2];    
     const dfloat rn3 = vr[sv3];    
     const dfloat sn3 = vs[sv3];    

     const dfloat sxe1 = -0.5*(rn1+sn1)*xe1 + 0.5*(1+rn1)*xe2 + 0.5*(1+sn1)*xe3;
     const dfloat sye1 = -0.5*(rn1+sn1)*ye1 + 0.5*(1+rn1)*ye2 + 0.5*(1+sn1)*ye3;

     const dfloat sxe2 = -0.5*(rn2+sn2)*xe1 + 0.5*(1+rn2)*xe2 + 0.5*(1+sn2)*xe3;
     const dfloat sye2 = -0.5*(rn2+sn2)*ye1 + 0.5*(1+rn2)*ye2 + 0.5*(1+sn2)*ye3;

     const dfloat sxe3 = -0.5*(rn1+sn3)*xe1 + 0.5*(1+rn3)*xe2 + 0.5*(1+sn3)*xe3;
     const dfloat sye3 = -0.5*(rn1+sn3)*ye1 + 0.5*(1+rn3)*ye2 + 0.5*(1+sn3)*ye3;
      // 1/ Area of triangle
     const dfloat vol = 0.5*((sxe2-sxe1)*(sye3-sye1) - (sxe3-sxe1)*(sye2-sye1));
      vgeo[Nvgeo*elm + IVID] = 1.0 / vol;
    }
  }
  // //Use mesh halo exchange 
  // mesh.halo->Exchange(vgeo, Nvgeo*Nsubcells, ogs_dfloat);

  
  // Compute Surface Geometric Factors
  sgeo = (dfloat*) calloc((mesh.Nelements+mesh.totalHaloPairs)*Nsubcells*Nfaces*Nsgeo, 
                           sizeof(dfloat));

  for(dlong e=0; e<(mesh.Nelements + mesh.totalHaloPairs); e++){
    dlong id = e*Nverts;
    const dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
    const dfloat xe2 = mesh.EX[id+1];
    const dfloat xe3 = mesh.EX[id+2];

    const dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
    const dfloat ye2 = mesh.EY[id+1];
    const dfloat ye3 = mesh.EY[id+2];

    for(int s= 0; s<Nsubcells; s++){
      const dlong elm = e*Nsubcells + s;  
      // Get all vertices
      const int sv1 = mEToV[s*Nverts + 0]; 
      const int sv2 = mEToV[s*Nverts + 1]; 
      const int sv3 = mEToV[s*Nverts + 2];

       const dfloat rn1 = vr[sv1];    
       const dfloat sn1 = vs[sv1]; 

       const dfloat rn2 = vr[sv2];    
       const dfloat sn2 = vs[sv2]; 

       const dfloat rn3 = vr[sv3];    
       const dfloat sn3 = vs[sv3];    

       const dfloat sxe1 = -0.5*(rn1+sn1)*xe1 + 0.5*(1+rn1)*xe2 + 0.5*(1+sn1)*xe3;
       const dfloat sye1 = -0.5*(rn1+sn1)*ye1 + 0.5*(1+rn1)*ye2 + 0.5*(1+sn1)*ye3;

       const dfloat sxe2 = -0.5*(rn2+sn2)*xe1 + 0.5*(1+rn2)*xe2 + 0.5*(1+sn2)*xe3;
       const dfloat sye2 = -0.5*(rn2+sn2)*ye1 + 0.5*(1+rn2)*ye2 + 0.5*(1+sn2)*ye3;

       const dfloat sxe3 = -0.5*(rn1+sn3)*xe1 + 0.5*(1+rn3)*xe2 + 0.5*(1+sn3)*xe3;
       const dfloat sye3 = -0.5*(rn1+sn3)*ye1 + 0.5*(1+rn3)*ye2 + 0.5*(1+sn3)*ye3;
  
        // face 1
        dfloat nx1 =  sye2-sye1;
        dfloat ny1 = -(sxe2-sxe1);
        dfloat  d1 = sqrt((nx1)*(nx1)+(ny1)*(ny1));

        sgeo[elm*Nfaces*Nsgeo + 0*Nsgeo + FXID] = 0.5*(sxe1 + sxe2);  
        sgeo[elm*Nfaces*Nsgeo + 0*Nsgeo + FYID] = 0.5*(sye1 + sye2);  
        sgeo[elm*Nfaces*Nsgeo + 0*Nsgeo + NXID] = nx1/d1;  
        sgeo[elm*Nfaces*Nsgeo + 0*Nsgeo + NYID] = ny1/d1;  
        sgeo[elm*Nfaces*Nsgeo + 0*Nsgeo + SAID] = d1; 

        // face2
        dfloat nx2 = sye3-sye2;
        dfloat ny2 = -(sxe3-sxe2);
        dfloat  d2 = sqrt((nx2)*(nx2)+(ny2)*(ny2));

        sgeo[elm*Nfaces*Nsgeo + 1*Nsgeo + FXID] = 0.5*(sxe2 + sxe3);  
        sgeo[elm*Nfaces*Nsgeo + 1*Nsgeo + FYID] = 0.5*(sye2 + sye3);  
        sgeo[elm*Nfaces*Nsgeo + 1*Nsgeo + NXID] = nx2/d2;  
        sgeo[elm*Nfaces*Nsgeo + 1*Nsgeo + NYID] = ny2/d2;  
        sgeo[elm*Nfaces*Nsgeo + 1*Nsgeo + SAID] = d2;  

        // face3
        dfloat nx3 =   sye1-sye3;
        dfloat ny3 = -(sxe1-sxe3);
        dfloat  d3 = sqrt((nx3)*(nx3)+(ny3)*(ny3));
        sgeo[elm*Nfaces*Nsgeo + 2*Nsgeo + FXID] = 0.5*(sxe1 + sxe3);  
        sgeo[elm*Nfaces*Nsgeo + 2*Nsgeo + FYID] = 0.5*(sye1 + sye3);  
        sgeo[elm*Nfaces*Nsgeo + 2*Nsgeo + NXID] = nx3/d3;  
        sgeo[elm*Nfaces*Nsgeo + 2*Nsgeo + NYID] = ny3/d3;  
        sgeo[elm*Nfaces*Nsgeo + 2*Nsgeo + SAID] = d3;  
        
    }
  }

  // Create Boundary Info: Might be needed
  // for(dlong e=0; e<(mesh.Nelements+mesh.totalHaloPairs); e++){
  for(dlong e=0; e<mesh.Nelements; e++){
    for(int f = 0; f<Nfaces; f++){
      int bc = mesh.EToB[e*mesh.Nfaces + f]; 

      if(bc){
        // get all elemnts and faces on this face
        for(int n=0; n<N; n++){
          // local element and face info 
          const int sem  =  mFToE[f*N+ n];
          const int sfm  =  mFToF[f*N+ n];
          const dlong id = e*Nsubcells*Nfaces + sem*Nfaces + sfm; 
          sgeo[id*Nsgeo + BCID] = bc; 
        }    
      }
    }
  }
}


void subcellTri2D::SetupOperators(){
// First Create Projection Matrix
int  cubN  = mesh.N+1; // Just to make sure it is dealiased
int  cubNp = 0; 
dfloat *cubr, *cubs, *cubw; 
mesh.CubatureNodesTri2D(cubN, &cubNp, &cubr, &cubs, &cubw); 

PM = (dfloat*)malloc(Nsubcells*mesh.Np*sizeof(dfloat));

dfloat *Ptemp = (dfloat*)malloc(Nsubcells*mesh.Np*sizeof(dfloat));

dfloat *_cx = (dfloat *)malloc(cubNp*sizeof(dfloat));  
dfloat *_cy = (dfloat *)malloc(cubNp*sizeof(dfloat));  

for(int s=0; s<Nsubcells; s++){
  
  // First compute cubature nodes on this element
  const int v1 = mEToV[s*Nverts + 0]; 
  const int v2 = mEToV[s*Nverts + 1]; 
  const int v3 = mEToV[s*Nverts + 2];

  const dfloat xe1 = vr[v1], ye1 = vs[v1];    
  const dfloat xe2 = vr[v2], ye2 = vs[v2];    
  const dfloat xe3 = vr[v3], ye3 = vs[v3];

  // Cubature nodes for this subcell
  for(int n=0; n<cubNp; n++){
    _cx[n] = -0.5*(cubr[n]+cubs[n])*xe1 + 0.5*(1+cubr[n])*xe2 + 0.5*(1+cubs[n])*xe3;
    _cy[n] = -0.5*(cubr[n]+cubs[n])*ye1 + 0.5*(1+cubr[n])*ye2 + 0.5*(1+cubs[n])*ye3;
  }

  int sk=0; 
  for(int i=0; i<mesh.N+1; i++){
      for(int j=0; j<mesh.N+1-i; j++){
        dfloat phi = 0;
        // integrate this basis function on the subcell s 
        for(int n=0; n<cubNp; n++){
          dfloat pn = 0; 
          mesh.OrthonormalBasisTri2D(_cx[n], _cy[n], i, j, &pn);
          phi += mJ[s]*cubw[n]*pn; 
        } 
        // Divide by the volume of subcell now
        Ptemp[s*mesh.Np + sk] =  1.0/(2.0*mJ[s])*phi; // J's are not needed for this affine mapping
        sk++;
    }
  }

}

// Get Vandermonde Matrix 
dfloat *V = (dfloat *) malloc(mesh.Np*mesh.Np*sizeof(dfloat));
mesh.VandermondeTri2D(mesh.N, mesh.Np, mesh.r, mesh.s, V);
// PM = Pt*V^-1
matrixRightSolve(Nsubcells, mesh.Np, Ptemp, mesh.Np, mesh.Np, V, PM); 
// Pseudo Inverse of PM to recontruct solution from FV
// This is a overdetermined system which needs LS approximation 
RM = (dfloat*)malloc(Nsubcells*mesh.Np*sizeof(dfloat));
for(int n=0; n<Nsubcells*mesh.Np; n++) 
 RM[n] = PM[n];

matrixPInverse(Nsubcells, mesh.Np, RM);

// Create conservative face interpolation and its inverse
// Required for DG-FV interface flux evaluation
int intNfp = cubN+1;
dfloat *intr = (dfloat *) malloc(intNfp*sizeof(dfloat));
dfloat *intw = (dfloat *) malloc(intNfp*sizeof(dfloat));
mesh.JacobiGQ(0, 0, cubN, intr, intw);

// 
dfloat *_fgr     = (dfloat *)malloc(intNfp*sizeof(dfloat));  
dfloat *ptemp0   = (dfloat *)malloc(N*mesh.Nfp*sizeof(dfloat));  
dfloat *V1D      = (dfloat *)malloc(mesh.Nfp*mesh.Nfp*sizeof(dfloat));
dfloat *r1D      = (dfloat *)malloc(mesh.Nfp*sizeof(dfloat));
dfloat *ptemp1   = (dfloat *)malloc(N*mesh.Nfp*sizeof(dfloat));  
dfloat *rtemp1   = (dfloat *)malloc(N*mesh.Nfp*sizeof(dfloat));  

// PFM = (dfloat *)malloc(N*mesh.Nfp*Nfaces*sizeof(dfloat));  
// RFM = (dfloat *)malloc(N*mesh.Nfp*Nfaces*sizeof(dfloat));  

PFM = (dfloat *)malloc(N*mesh.Nfp*sizeof(dfloat));  
RFM = (dfloat *)malloc(N*mesh.Nfp*sizeof(dfloat));  

// for(int f=0; f<Nfaces; f++){
for(int f=0; f<1; f++){
// Create DG operator
 dfloat *rFace;
 if (f==0) rFace = mesh.r;
 if (f==1) rFace = mesh.r;
 if (f==2) rFace = mesh.s;

  for (int i=0;i<mesh.Nfp;i++)
    r1D[i] = rFace[mesh.faceNodes[f*mesh.Nfp+i]];

 mesh.Vandermonde1D(mesh.N, mesh.Nfp, r1D, V1D); 

for(int s=0; s<N; s++){
  // local element and face info 
  const int sem  =  mFToE[f*N+s];
  const int sfm  =  mFToF[f*N+s];
// printf("%d %d \n", sem, sem*Nverts + faceVertices[sfm*NfaceVertices+0]);
  const int vid1 = mEToV[sem*Nverts + faceVertices[sfm*NfaceVertices+0]];
  const int vid2 = mEToV[sem*Nverts + faceVertices[sfm*NfaceVertices+1]];
  //
// printf("I am here\n");
  dfloat rn1 = 0, rn2 = 0;
  if(f==0){rn1 = vr[vid1];  rn2 = vr[vid2];} 
  if(f==1){rn1 = vr[vid1];  rn2 = vr[vid2];} 
  if(f==2){rn1 = vs[vid1];  rn2 = vs[vid2];} 
  //
  // A_subcellFace/A_referenceFace
  const dfloat jac = (rn2-rn1)/2.0;

   for(int n=0; n<intNfp; n++){
    _fgr[n] = rn1 + 0.5*(1.0 + intr[n])*(rn2 - rn1);
   }

  //
  int sk=0; 
  for(int n=0; n<mesh.Nfp; n++){
    dfloat phi = 0; 
    // printf("f = %d id = %d\n", f, mesh.faceNodes[f*mesh.Nfp + n]);
    for(int i=0; i<intNfp; i++){
      dfloat pn = 0;
      mesh.OrthonormalBasis1D(_fgr[i], n, &pn);
      phi += jac*intw[i]*pn; 
    }
    ptemp0[s*mesh.Nfp + sk] =  1.0/(2.0*jac)*phi; 
    sk++;
   }
  }
 
  matrixRightSolve(N, mesh.Nfp,ptemp0,mesh.Nfp, mesh.Nfp, V1D, ptemp1);

  for(int n=0; n<N*mesh.Nfp; n++)
    rtemp1[n] = ptemp1[n]; 

  matrixPInverse(N, mesh.Nfp, rtemp1);

  // fill operators
  for(int i=0; i<N; i++){
    for(int j=0; j<mesh.Nfp; j++){
      PFM[f*mesh.Nfp*N + i*mesh.Nfp + j] = ptemp1[i*mesh.Nfp + j];  
      RFM[f*mesh.Nfp*N + j*N + i]        = rtemp1[j*N + i];  
    }
  }
}

// Compute Subcell Lift Matrix

SLIFT = (dfloat *)malloc(mesh.Np*Nfaces*N*sizeof(dfloat)); 

for(int f=0; f<Nfaces; f++){
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<N; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Nfp; m++){
        const int id = f*mesh.Nfp + m; 
        sum += mesh.LIFT[i*mesh.Nfp*mesh.Nfaces + id]*RFM[m*N + j]; 
      }
     SLIFT[i*(Nfaces*N) + (f*N + j)] = sum; 
    }
  }
}


  // for(int i=0; i<mesh.Np; i++){
  //   for(int j=0; j<(N*Nfaces); j++){
  //    printf("%.4e ",SLIFT[i*(N*Nfaces) + j] );
  //   }
  //   printf("\n");
  // }
 
free(_fgr); free(r1D);      free(V1D); 
free(ptemp0); free(ptemp1); free(rtemp1); 
free(_cx);  free(_cy); 
free(cubr); free(cubs); free(cubw);
free(V);


}
