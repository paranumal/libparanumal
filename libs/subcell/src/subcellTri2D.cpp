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

subcellTri2D::subcellTri2D(mesh_t& _mesh, settings_t& _settings):
   subcell_t(_mesh, _settings) {
   _settings.getSetting("SUBCELL NUMBER", N);
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
  // ElementList = (dfloat *) calloc(mesh.Nelements, sizeof(dfloat)); 
  ElementList = (dlong *) calloc(mesh.Nelements, sizeof(dlong)); 


  // N = N_;
  // Nfp = N+1;
  // Np = (N+1)*(N+2)/2;

  // /* Nodal Data */
  // r = (dfloat *) malloc(Np*sizeof(dfloat));
  // s = (dfloat *) malloc(Np*sizeof(dfloat));
  // NodesTri2D(N, r, s);

  // faceNodes = (int *) malloc(Nfaces*Nfp*sizeof(int));
  // FaceNodesTri2D(N, r, s, faceNodes);

  // vertexNodes = (int*) calloc(Nverts, sizeof(int));
  // VertexNodesTri2D(N, r, s, vertexNodes);

  // dfloat *V = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // VandermondeTri2D(N, Np, r, s, V);

  // //Mass matrix
  // MM = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // MassMatrixTri2D(Np, V, MM);
  // free(V);

  // //D matrices
  // Dr = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // Ds = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // DmatrixTri2D(N, Np, r, s, Dr, Ds);

  // // Weak D matrices
  // DWr = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // DWs = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  // DWmatrixTri2D(N, Np, r, s, MM, DWr, DWs);

  // LIFT = (dfloat *) malloc(Np*Nfaces*Nfp*sizeof(dfloat));
  // LIFTmatrixTri2D(N, faceNodes, r, s, LIFT);

  // /* Plotting data */
  // int plotN = N + 3; //enriched interpolation space for plotting
  // plotNp = (plotN+1)*(plotN+2)/2;

  // /* Plotting nodes */
  // plotR = (dfloat *) malloc(plotNp*sizeof(dfloat));
  // plotS = (dfloat *) malloc(plotNp*sizeof(dfloat));
  // EquispacedNodesTri2D(plotN, plotR, plotS);

  // plotNelements = plotN*plotN;
  // plotNverts = 3;
  // plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));
  // EquispacedEToVTri2D(plotN, plotEToV);

  // plotInterp = (dfloat *) malloc(Np*plotNp*sizeof(dfloat));
  // InterpolationMatrixTri2D(N, Np, r, s, plotNp, plotR, plotS, plotInterp);
}


void subcellTri2D::OccaSetup(){
 occa::properties kernelInfo = props; //copy base properties
 
 // Element List
 o_ElementList = device.malloc(mesh.Nelements*sizeof(dlong), ElementList);
 // o_ElementList = device.malloc(mesh.Nelements*sizeof(dfloat), ElementList);
 o_invVT       = device.malloc(mesh.Np*mesh.Np*sizeof(dfloat), invVT); 
 o_ModMap      = device.malloc(mesh.Np*sizeof(int), ModeMap); 
 o_LSF         = device.malloc(mesh.N*sizeof(dfloat), LSF); 


 skylineKernel = buildKernel(device, SUBCELL_DIR "/okl/"
                                        "detectorTri2D.okl",
                                        "skylineTri2D",
                                        kernelInfo, comm);

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
  // xf = (dfloat *) malloc(mesh.Nelements*Nsubcells*Nfaces*sizeof(dfloat));
  // yf = (dfloat *) malloc(mesh.Nelements*Nsubcells*Nfaces*sizeof(dfloat));

  // xe = (dfloat *) malloc(mesh.Nelements*Nsubcells*sizeof(dfloat));
  // ye = (dfloat *) malloc(mesh.Nelements*Nsubcells*sizeof(dfloat));

  // dlong cnt = 0; 
  // for(dlong e=0; e<mesh.Nelements; e++){
  //   dlong id = e*Nverts;
  //   dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
  //   dfloat xe2 = mesh.EX[id+1];
  //   dfloat xe3 = mesh.EX[id+2];

  //   dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
  //   dfloat ye2 = mesh.EY[id+1];
  //   dfloat ye3 = mesh.EY[id+2];

  //   for(int s= 0; s<Nsubcells; s++){
  //     xe[e*Nsubcells+s] = (xe1 + xe2 + xe3)/3.0;
  //     ye[e*Nsubcells+s] = (ye1 + ye2 + ye3)/3.0;
  //     //
  //     for(int f =0; f<Nfaces; f++){
  //       dfloat rn = fr[s*Nfaces+f]; 
  //       dfloat sn = fs[s*Nfaces+f]; 

  //       xf[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
  //       yf[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
  //       ++cnt;
  //     }
  //   }
  // }




}
