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

#include "stab.hpp"

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

void stab_t::stabSetupSubcellTri2D(){

// Create minor grid i.e. EToV, (cr,cs), Jacobian etc.
CellCreateMinorGridTri2D(); 

// Create Local Connection i.e. reference sub-tesselation
CellLocalConnectTri2D(); 

// Form global connectivities on subcell level
CellGlobalConnectTri2D(); 

// Compute Geometric Factors
CellGeometricFactorsTri2D(); 

// Compute Geometric Factors
CellSetupOperatorsTri2D(); 

// // printf("Here Computed GeometricFactor on Subcells\n");

// int blockMax = 256;
// if (platform.device.mode() == "CUDA") blockMax = 512;

// props["defines/" "s_N"]         = N;
// props["defines/" "s_Np"]        = Np;
// props["defines/" "s_Ncells"]    = Nsubcells;
// props["defines/" "s_Nfaces"]    = Nfaces;
// props["defines/" "s_NfacesNfp"] = Nfaces*N;   


// props["defines/" "s_Nvgeo"]= Nvgeo;
// props["defines/" "s_Nsgeo"]= Nsgeo;

// props["defines/" "s_CXID"]= CXID;
// props["defines/" "s_CYID"]= CYID;
// props["defines/" "s_IVID"]= IVID;

// props["defines/" "s_FXID"]= FXID;
// props["defines/" "s_FYID"]= FYID;
// props["defines/" "s_NXID"]= NXID;
// props["defines/" "s_NYID"]= NYID;
// props["defines/" "s_SAID"]= SAID;
// props["defines/" "s_BCID"]= BCID;

// int maxNodesV = std::max(mesh.Np, Nsubcells);
// props["defines/" "s_maxNodesV"]  = maxNodesV;

// int maxNodesS = std::max(Nfaces*mesh.Nfp, Nfaces*N);
// props["defines/" "s_maxNodesS"]  = maxNodesS;


// // int NblockV = std::max(1, blockMax/maxNodesV);
// int NblockV = 1; 
// props["defines/" "s_NblocksV"]= NblockV;

// int NblockS = 1;
// // int NblockS = std::max(1, blockMax/maxNodesS);
// props["defines/" "s_NblocksS"]= NblockS;


// std::string oklFilePrefix = STAB_DIR "/okl/";
// std::string oklFileSuffix = ".okl";

// std::string fileName, kernelName;
// fileName      = oklFilePrefix + "subcell" + oklFileSuffix;
// kernelName    = "projectFVTri2D";
// projectFVKernel  = platform.buildKernel(fileName, kernelName, props);

// // kernelName    = "projectViscosity" + suffix;
// kernelName    = "projectViscosity";
// projectViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

// // kernelName    = "projectViscosity" + suffix;
// kernelName    = "maskElements";
// maskElementsKernel  = platform.buildKernel(fileName, kernelName, props);





}


void stab_t::CellFindBestMatchTri2D(dfloat x1, dfloat y1, dlong eP, int fP, 
                                    memory<int> &elist,  memory<dfloat> &x2,  memory<dfloat> &y2, 
                                    int &nE, int &nF){

  int eIndex = 0;
  int fIndex = 0;
  dfloat mindist2 = 1e12;  

  const dlong shift = eP*Nsubcells*Nfaces; 

  for(int e=0;e<N;++e){
    /* next element */
    const int e2 = elist[e + fP*N];
    /* check faces */
    for(int f2=0; f2<Nfaces; f2++){  
      /* distance between target and next node */
      const dfloat dist2 = pow(x1-x2[e2*Nfaces + f2 + shift],2) + 
                           pow(y1-y2[e2*Nfaces + f2 + shift],2);
          /* if next node is closer to target update match */
      if(dist2<mindist2){
        mindist2 = dist2;
        eIndex = e2; 
        fIndex = f2;
      }

    }
  }
   
  if(mindist2>1e-6) { LIBP_FORCE_ABORT("bad cell connection ");}
  
  nE = eIndex; 
  nF = fIndex; 
}


// structure used to encode vertices that make
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct {
  hlong v[4]; // vertices on face
  dlong element, elementN;
  int face, faceN;    // face info
  int rank, rankN; // N for neighbor face info

}face_t;


// Create Local Connectivity Inside 
void stab_t::CellLocalConnectTri2D(){

  // local cell to subcell trace info   
  mEToE.malloc(Nsubcells*Nfaces,0);
  mEToF.malloc(Nsubcells*Nfaces,0);
  /* build list of faces */
  memory<face_t> faces(Nsubcells*Nfaces);

  // #pragma omp parallel for collapse(2)
  for(dlong e=0;e<Nsubcells;++e){
    for(int f=0;f<Nfaces;++f){
      const dlong id = f + e*Nfaces;

      for(int n=0;n<NfaceVertices;++n){
        dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
        faces[id].v[n] = mEToV[vid];
      }

      std::sort(faces[id].v, faces[id].v+NfaceVertices, std::less<hlong>());

      faces[id].element = e;
      faces[id].face = f;

      faces[id].elementN= -1;
      faces[id].faceN = -1;
    }
  }

  /* sort faces by their vertex number pairs */
  sort(faces.ptr(), faces.ptr()+Nsubcells*Nfaces,
       [&](const face_t& a, const face_t& b) {
         return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                             b.v, b.v+NfaceVertices);
       });

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  // #pragma omp parallel for
  for(dlong cnt=0;cnt<Nsubcells*Nfaces-1;++cnt){
    if(std::equal(faces[cnt].v, faces[cnt].v+NfaceVertices,
                  faces[cnt+1].v)){
      // match
      faces[cnt].elementN = faces[cnt+1].element;
      faces[cnt].faceN    = faces[cnt+1].face;

      faces[cnt+1].elementN = faces[cnt].element;
      faces[cnt+1].faceN    = faces[cnt].face;
    }
  }


   /* resort faces back to the original element/face ordering */
  sort(faces.ptr(), faces.ptr()+Nsubcells*Nfaces,
       [](const face_t& a, const face_t& b) {
         if(a.element < b.element) return true;
         if(a.element > b.element) return false;

         return (a.face < b.face);
       });


  /* extract the element to element and element to face connectivity */
  Nint =0, Next=0; 
  // #pragma omp parallel for collapse(2)
  for(dlong e=0;e<Nsubcells;++e){
    int eext=0; 
    for(int f=0;f < Nfaces;++f){
      const dlong id = f + e*Nfaces;
      mEToE[id] = faces[id].elementN;
      mEToF[id] = faces[id].faceN;
      if( mEToE[e*Nfaces + f]< 0 ){eext = 1;} 
    }

    if(eext){      
      Next++;
    }else{
      Nint++;
    }
  }

  faces.free();

  ielist.malloc(Nint,0);
  eelist.malloc(Next,0);
  // 
  int ske=0, ski=0; 
  for(dlong e=0;e<Nsubcells;++e){
    int eext = 0; 
    for(int f=0;f<Nfaces;++f){
      if(mEToE[e*Nfaces + f]<0){eext = 1;}
    }

    if(eext)
      eelist[ske++] = e;        
    else
      ielist[ski++] = e; 
  }

// external connection 
  dfloat deps = 1.;
  while((1.+deps)>1.) deps *= 0.5;
  const dfloat NODETOL = 1000.*deps;
  
  memory<int> fcnt(Nfaces,0); 
  // local cell to subcell trace info   
  mFToE.malloc(N*Nfaces,0);
  mFToF.malloc(N*Nfaces,0);

  for(dlong e=0;e<Nsubcells;++e){
    for(int f=0;f<Nfaces;++f){
      if(mEToE[e*Nfaces + f]<0){ // local boundary element

        dfloat rf = 0, sf = 0; 
        for(int n=0;n<NfaceVertices;++n){
          dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
          rf += vr[mEToV[vid]]/NfaceVertices; 
          sf += vs[mEToV[vid]]/NfaceVertices;
        }

        if(mesh.elementType==Mesh::TRIANGLES){
          if(fabs(sf+1)<NODETOL){
            mFToE[0*N+ fcnt[0]] = e;
            mFToF[0*N+ fcnt[0]] = f;
            fcnt[0]++; 
          }
          if(fabs(rf+sf)<NODETOL){
            mFToE[1*N+ fcnt[1]] = e;
            mFToF[1*N+ fcnt[1]] = f;
            fcnt[1]++;
          }
          if(fabs(rf+1)<NODETOL){
            mFToE[2*N+fcnt[2]] = e;
            mFToF[2*N+fcnt[2]] = f;
            fcnt[2]++;
          }
        }
      }
    }
  }

o_mFToE = platform.malloc<int>(mFToE); 
o_mFToF = platform.malloc<int>(mFToF); 

o_eelist  = platform.malloc<int>(eelist); 
o_ielist  = platform.malloc<int>(ielist); 
}


void stab_t::CellGlobalConnectTri2D(){

  const dlong Nelements = mesh.Nelements; 
  emapP.malloc(mesh.Nelements*Nsubcells*Nfaces,0);
  fmapP.malloc(mesh.Nelements*Nsubcells*Nfaces,0);

  // first connect internal elements
  for(dlong e=0; e<Nelements; e++){
    // for all cells
    for(int s =0; s<Nsubcells; s++){
      const dlong eshift = e*Nsubcells*Nfaces + s*Nfaces; 
      for(int f=0; f<Nfaces; f++){
        // check local connectivity
        int ep = mEToE[s*Nfaces + f]; 
        int fp = mEToF[s*Nfaces + f]; 
        if(!(ep<0)){ // local connection, we dont need to hold will check later!!!!
          emapP[eshift + f] = ep + e*Nsubcells; 
          fmapP[eshift + f] = fp; //  
        }
      }
    }
  }

  // //check if we're using a periodic box mesh
  // int periodicFlag = 0;
  // if (mesh.settings.compareSetting("MESH FILE","BOX") &&
  //     mesh.settings.compareSetting("BOX BOUNDARY FLAG","-1"))
  //   periodicFlag = 1;

  // //box dimensions
  // dfloat DIMX, DIMY;
  // mesh.settings.getSetting("BOX DIMX", DIMX);
  // mesh.settings.getSetting("BOX DIMY", DIMY);

  // //box is centered at the origin
  // DIMX /= 2.0;
  // DIMY /= 2.0;

  // Compute Face Centers and Connect
  memory<dfloat> xf((mesh.Nelements + mesh.totalHaloPairs)*Nsubcells*Nfaces,0.0);
  memory<dfloat> yf((mesh.Nelements + mesh.totalHaloPairs)*Nsubcells*Nfaces,0.0);

  dlong cnt = 0; 
  for(dlong e=0; e<mesh.Nelements; e++){
  // for(dlong e=0; e<(mesh.Nelements + mesh.totalHaloPairs); e++){
    dlong id = e*Nverts;
    dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh.EX[id+1];
    dfloat xe3 = mesh.EX[id+2];

    dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh.EY[id+1];
    dfloat ye3 = mesh.EY[id+2];

    for(int s= 0; s<Nsubcells; s++){
      //
      for(int f =0; f<Nfaces; f++){
        dfloat rn = fr[s*Nfaces+f]; 
        dfloat sn = fs[s*Nfaces+f]; 

        xf[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
        yf[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
        ++cnt;
      }
    }
  }

  mesh.halo.Exchange(xf, Nsubcells*Nfaces);
  mesh.halo.Exchange(yf, Nsubcells*Nfaces);

  for(dlong e=0; e<Nelements; e++){
    for(int f=0;f<Nfaces;++f){
      dlong eP = mesh.EToE[e*Nfaces + f]; 
      int fP   = mesh.EToF[e*Nfaces + f]; 

      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e; fP = f;
      }

      dfloat offsetX = 0.0;
      dfloat offsetY = 0.0;

      // if (periodicFlag) {
      //   //if the mesh is periodic, this is more complicated.
      //   // check if this face is on a boundary face
      //   bool top=true, bottom=true, left=true, right=true;
      //   for(int n=0;n<NfaceVertices;++n){
      //     dlong vid = e*mesh.Nverts + mesh.faceVertices[f*NfaceVertices+n];
      //     if (fabs(mesh.EX[vid]-DIMX)>1e-4) right = false;
      //     if (fabs(mesh.EX[vid]+DIMX)>1e-4) left = false;
      //     if (fabs(mesh.EY[vid]-DIMY)>1e-4) top = false;
      //     if (fabs(mesh.EY[vid]+DIMY)>1e-4) bottom = false;
      //   }
      //   if (right)  offsetX = -2.0*DIMX;
      //   if (left)   offsetX =  2.0*DIMX;
      //   if (top)    offsetY = -2.0*DIMY;
      //   if (bottom) offsetY =  2.0*DIMY;
      // }

      // for each subcell at this face find neighboor element
      for(int n=0; n<N; n++){
        // local element and face info 
        const int sem =  mFToE[f*N+ n];
        const int sfm =  mFToF[f*N+ n];
        //
        const dlong idM = e*Nsubcells*Nfaces + sem*Nfaces + sfm; 
        dfloat xM = xf[idM] + offsetX;
        dfloat yM = yf[idM] + offsetY;


        // printf("%d %d %d %d \n", mesh.rank, e, sem, sfm);
        int idE, idF; // mEToE[sem*Nfaces + sfm]; 
        CellFindBestMatchTri2D(xM, yM, eP, fP, mFToE, xf, yf, idE, idF);

        const dlong eshift = e*Nsubcells*Nfaces + sem*Nfaces + sfm;  
        emapP[eshift]     = idE + eP*Nsubcells; 
        fmapP[eshift]     = idF; //

      }
    }
  }

  o_emapP = platform.malloc<dlong>(emapP); 
  o_fmapP = platform.malloc<dlong>(fmapP); 
  
  // mesh.halo.Exchange(o_emapP, Nsubcells*Nfaces);
  // mesh.halo.Exchange(o_fmapP, Nsubcells*Nfaces);

//   if(mesh.rank==0){ 
//   for(dlong e=0; e<Nelements; e++){
//     for(int f=0;f<Nfaces;++f){
//       printf("%d  ", mesh.EToE[e*mesh.Nfaces + f]);
//       }
//       printf("\n");
//     }
//   }

 // printf("\n");
 //  printf("\n");
 // if(mesh.rank==1){ 
 //  for(dlong e=0; e<Nelements; e++){
 //    for(int f=0;f<Nfaces;++f){
 //      printf("%d  ", mesh.EToE[e*mesh.Nfaces + f]);
 //      }
 //      printf("\n");
 //    }
 //  }


 //  printf("\n");
 //  printf("\n");
 // if(mesh.rank==1){ 
 //  for(dlong e=0; e<Nelements; e++){
 //    for(int s=0; s<Nsubcells; s++){
 //        for(int f=0;f<Nfaces;++f){
 //      printf("%d  ", emapP[e*Nsubcells*mesh.Nfaces + s*mesh.Nfaces + f]);
 //      }
 //      printf("\n");
 //    }
 //    printf("\n");
 //  }
 //  }


  // if(mesh.rank==0){ 
  // for(dlong e=0; e<Nelements; e++){
  //   for(int f=0;f<Nfaces;++f){
  //     printf("%d  ", mesh.EToV[e*mesh.Nverts + f]);
  //     }
  //     printf("\n");
  //   }
  // }


// printf("\n");
//   printf("\n");

  // if(mesh.rank==0){ 
  // for(dlong e=0; e<Nelements; e++){
  //   for(int f=0;f<Nfaces;++f){
  //      const dlong id = e*mesh.Nfaces*mesh.Nfp + f*mesh.Nfp + 0; 
  //       const dlong vid = mesh.vmapP[id]; 
  //       const dlong ep  = vid/mesh.Np;     
  //     printf("%d  ", ep);
  //     }
  //     printf("\n");
  //   }
  // }

//  printf("\n");
//   printf("\n");
//  if(mesh.rank==1){ 
//   for(dlong e=0; e<Nelements; e++){
//     for(int f=0;f<Nfaces;++f){
//        const dlong id = e*mesh.Nfaces*mesh.Nfp + f*mesh.Nfp + 0; 
//         const dlong vid = mesh.vmapP[id]; 
//         const dlong ep  = vid/mesh.Np;     
//       printf("%d  ", ep);
//       }
//       printf("\n");
//     }
//   }

}


void stab_t::CellGeometricFactorsTri2D(){
  Nvgeo = 3; 
  Nsgeo = 6; 
  //
  CXID  = 0;
  CYID  = 1;
  IVID  = 2;

  FXID  = 0;
  FYID  = 1;
  NXID  = 2;
  NYID  = 3;
  SAID  = 4;
  BCID  = 5;

  vgeo.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nsubcells*Nvgeo);  

  // Compute Volume Geometric Facors
  // for(dlong e=0; e<(mesh.Nelements+mesh.totalHaloPairs); e++){
  for(dlong e=0; e<mesh.Nelements; e++){
    dlong id = e*Nverts;
    const dfloat xe1 = mesh.EX[id+0]; /* x-coordinates of vertices */
    const dfloat xe2 = mesh.EX[id+1];
    const dfloat xe3 = mesh.EX[id+2];

    const dfloat ye1 = mesh.EY[id+0]; /* y-coordinates of vertices */
    const dfloat ye2 = mesh.EY[id+1];
    const dfloat ye3 = mesh.EY[id+2];

    for(int s =0; s<Nsubcells; s++){
      // local cell id 
      const dlong elm = e*Nsubcells + s;  
      // Compute Cell Centers
      const dfloat rn = cr[s]; 
      const dfloat sn = cs[s]; 
      // Cell centers
      vgeo[Nvgeo*elm + CXID] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      vgeo[Nvgeo*elm + CYID] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;

      // Get all vertices
      const int sv1 = mEToV[s*Nverts + 0]; 
      const int sv2 = mEToV[s*Nverts + 1]; 
      const int sv3 = mEToV[s*Nverts + 2];

      const dfloat rn1 = vr[sv1],sn1 = vs[sv1];    
      const dfloat rn2 = vr[sv2],sn2 = vs[sv2];    
      const dfloat rn3 = vr[sv3],sn3 = vs[sv3];    

      const dfloat sxe1 = -0.5*(rn1+sn1)*xe1 + 0.5*(1+rn1)*xe2 + 0.5*(1+sn1)*xe3;
      const dfloat sye1 = -0.5*(rn1+sn1)*ye1 + 0.5*(1+rn1)*ye2 + 0.5*(1+sn1)*ye3;

      const dfloat sxe2 = -0.5*(rn2+sn2)*xe1 + 0.5*(1+rn2)*xe2 + 0.5*(1+sn2)*xe3;
      const dfloat sye2 = -0.5*(rn2+sn2)*ye1 + 0.5*(1+rn2)*ye2 + 0.5*(1+sn2)*ye3;

      const dfloat sxe3 = -0.5*(rn3+sn3)*xe1 + 0.5*(1+rn3)*xe2 + 0.5*(1+sn3)*xe3;
      const dfloat sye3 = -0.5*(rn3+sn3)*ye1 + 0.5*(1+rn3)*ye2 + 0.5*(1+sn3)*ye3;
      // 1/ Area of triangle
      const dfloat vol = 0.5*((sxe2-sxe1)*(sye3-sye1) - (sxe3-sxe1)*(sye2-sye1));
      vgeo[Nvgeo*elm + IVID] = 1.0 / vol;
    }
  }

 
  // Compute Surface Geometric Factors
  sgeo.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nsubcells*Nfaces*Nsgeo,0.0);

  // for(dlong e=0; e<(mesh.Nelements + mesh.totalHaloPairs); e++){
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

       const dfloat sxe3 = -0.5*(rn3+sn3)*xe1 + 0.5*(1+rn3)*xe2 + 0.5*(1+sn3)*xe3;
       const dfloat sye3 = -0.5*(rn3+sn3)*ye1 + 0.5*(1+rn3)*ye2 + 0.5*(1+sn3)*ye3;
  
        // face 1
        dfloat nx1 =   sye2-sye1; 
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

 // Get halo info 
 mesh.halo.Exchange(vgeo, Nvgeo*Nsubcells);
 mesh.halo.Exchange(sgeo, Nsgeo*Nsubcells*Nfaces);
 
 o_vgeo = platform.malloc<dfloat>(vgeo);
 o_sgeo = platform.malloc<dfloat>(sgeo);
}


void stab_t::CellSetupOperatorsTri2D(){
// First Create Projection Matrix
int  cubN  = mesh.N+1; // Just to make sure it is dealiased
int  cubNp = 0; 
memory<dfloat> cubr, cubs, cubw;  
mesh.CubatureNodesTri2D(cubN, cubNp, cubr, cubs, cubw);

memory<dfloat>_cx(cubNp, 0.0);  
memory<dfloat>_cy(cubNp, 0.0);  

memory<dfloat> Ptemp(Nsubcells*mesh.Np, 0.0);
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
  for(int i=0; i<(mesh.N+1); i++){
      for(int j=0; j<(mesh.N+1-i); j++){
        dfloat phi = 0;
        // integrate this basis function on the subcell s 
        for(int n=0; n<cubNp; n++){
          dfloat pn = 0; 
          mesh.OrthonormalBasisTri2D(_cx[n], _cy[n], i, j, pn);
          phi += mJ[s]*cubw[n]*pn; 
        } 
        // Divide by the volume of subcell now
        Ptemp[s*mesh.Np + sk] =  1.0/(2.0*mJ[s])*phi; // J's are not needed for this affine mapping
        sk++;
    }
  }
}

// Get Vandermonde Matrix 
memory<dfloat> V2D; 
mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, V2D);

// PM = Ptemp*V^-1
PM.malloc(Nsubcells*mesh.Np);
linAlg_t::matrixRightSolve(Nsubcells, mesh.Np, Ptemp, mesh.Np, mesh.Np, V2D, PM); 


RM.malloc(Nsubcells*mesh.Np);  RM.copyFrom(PM, Nsubcells*mesh.Np);  

/* Do not Forget!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
linAlg_t::matrixPseudoInverse(Nsubcells, mesh.Np, RM);
*/

// Create conservative face interpolation and its inverse
// Required for DG-FV interface flux evaluation
int intNfp = cubN+1;
memory<dfloat> intr, intw; 
mesh.JacobiGQ(0, 0, cubN, intr, intw);

memory<dfloat> _fgr  (intNfp);
// memory<dfloat> V1D   (mesh.Nfp*mesh.Nfp);
memory<dfloat> V1D;
memory<dfloat> r1D   (mesh.Nfp);

memory<dfloat> ptemp0(N*mesh.Nfp);
memory<dfloat> ptemp1(N*mesh.Nfp);
memory<dfloat> rtemp1(N*mesh.Nfp);

PFM.malloc(N*mesh.Nfp,0.0);
RFM.malloc(N*mesh.Nfp,0.0);

// for(int f=0; f<Nfaces; f++){
for(int f=0; f<1; f++){
// Create DG operator
 memory<dfloat> rFace(mesh.Np);
 if (f==0) rFace.copyFrom(mesh.r, mesh.Np);
 if (f==1) rFace.copyFrom(mesh.r, mesh.Np);
 if (f==2) rFace.copyFrom(mesh.s, mesh.Np);

  for (int i=0;i<mesh.Nfp;i++)
    r1D[i] = rFace[mesh.faceNodes[f*mesh.Nfp+i]];

 mesh.Vandermonde1D(mesh.N, r1D, V1D); 

for(int s=0; s<N; s++){
  // local element and face info 
  const int sem  =  mFToE[f*N+s];
  const int sfm  =  mFToF[f*N+s];

  const int vid1 = mEToV[sem*Nverts + faceVertices[sfm*NfaceVertices+0]];
  const int vid2 = mEToV[sem*Nverts + faceVertices[sfm*NfaceVertices+1]];

  dfloat rn1 = 0, rn2 = 0;
  if(f==0){rn1 = vr[vid1];  rn2 = vr[vid2];} 
  if(f==1){rn1 = vr[vid1];  rn2 = vr[vid2];} 
  if(f==2){rn1 = vs[vid1];  rn2 = vs[vid2];} 

  // // A_subcellFace/A_referenceFace
  const dfloat jac = (rn2-rn1)/2.0;

   for(int n=0; n<intNfp; n++){
    // printf("%.4f %.4f %.4f", intr[n], rn2, rn1);
    _fgr[n] = rn1 + 0.5*(1.0 + intr[n])*(rn2 - rn1);
   }

  //
  int sk=0; 
  for(int n=0; n<mesh.Nfp; n++){
    dfloat phi = 0; 
    for(int i=0; i<intNfp; i++){
      dfloat pn = 0;
      mesh.OrthonormalBasis1D(_fgr[i], n, pn);
      phi += jac*intw[i]*pn; 
    }
    ptemp0[s*mesh.Nfp + sk] =  1.0/(2.0*jac)*phi; 
    sk++;
   }
  }
 

  linAlg_t::matrixRightSolve(N, mesh.Nfp,ptemp0,mesh.Nfp, mesh.Nfp, V1D, ptemp1);
  
  rtemp1.copyFrom(ptemp1, N*mesh.Nfp);

  /* Do not forget ****************************************************************************
  linAlg_t::matrixPseudoInverse(N, mesh.Nfp, rtemp1);
  */

  // fill operators
  for(int i=0; i<N; i++){
    for(int j=0; j<mesh.Nfp; j++){
      PFM[f*mesh.Nfp*N + i*mesh.Nfp + j] = ptemp1[i*mesh.Nfp + j];  
      RFM[f*mesh.Nfp*N + j*N + i]        = rtemp1[j*N + i];  
    }
  }
}




// // Compute Subcell Lift Matrix
// SLIFT.malloc(mesh.Np*Nfaces*N);
// for(int f=0; f<Nfaces; f++){
//   for(int i=0; i<mesh.Np; i++){
//     for(int j=0; j<N; j++){
//       dfloat sum = 0; 
//       for(int m=0; m<mesh.Nfp; m++){
//         const int id = f*mesh.Nfp + m; 
//         sum += mesh.LIFT[i*mesh.Nfp*mesh.Nfaces + id]*RFM[m*N + j]; 
//       }
//      SLIFT[i*(Nfaces*N) + (f*N + j)] = sum; 
//     }
//   }
// }

// // Project Noddal Solution to Cell Vertices Directly
// memory<dfloat> PVMT = PVM;
// linAlg_t::matrixTranspose(Np, mesh.Np, PVM, mesh.Np, PVMT, Np); 
// o_PVM = platform.malloc<dfloat>(PVMT); 

// Project Nodal Solution to Cell Averages
memory<dfloat> PMT(mesh.Np*Nsubcells);
linAlg_t::matrixTranspose(Nsubcells, mesh.Np, PM, mesh.Np, PMT, Nsubcells); 
o_PM = platform.malloc<dfloat>(PMT); 

// Recontruct Nodal Solution From Cell Averages
memory<dfloat> RMT(mesh.Np*Nsubcells); 
linAlg_t::matrixTranspose(mesh.Np, Nsubcells, RM, Nsubcells, RMT, mesh.Np); 
o_RM = platform.malloc<dfloat>(RMT); 



/*
// Project Nodal Face Solution To Cell Face Averages
memory<dfloat> PFMT(N*mesh.Nfp);
linAlg_t::matrixTranspose(N, mesh.Nfp, PFM, mesh.Nfp, PFMT, N); 
o_PFM = platform.malloc<dfloat>(PFMT); 

// Recontruct Nodal Solution From Cell Face Averages
memory<dfloat> RFMT(mesh.Nfp*N);
linAlg_t::matrixTranspose(mesh.Nfp, N, RFM, N, RFMT, mesh.Nfp); 
o_RFM = platform.malloc<dfloat>(RFMT); 

*/





// // Lift Cell Averages Directly to Nodal Points 
// memory<dfloat> SLIFTT = SLIFT;
// linAlg_t::matrixTranspose(mesh.Np, (N*Nfaces), SLIFT, (N*Nfaces), SLIFTT, mesh.Np); 
// o_SLIFT = platform.malloc<dfloat>(SLIFTT); 

}



void stab_t::CellCreateMinorGridTri2D(){

  // Mesh face node ids
  mFaceNodes.malloc(mesh.Nfaces*mesh.Nfp,0); 
  mFaceNodes.copyFrom(mesh.faceNodes); 
  o_mFaceNodes = platform.malloc<int>(mFaceNodes); 

  // Using triangle does not have to be!!!!
  // Subcell could be in different topology but not now!
  Nverts = mesh.Nverts; 
  Nfaces = mesh.Nfaces; 
  // 
  NfaceVertices = mesh.NfaceVertices; 
  faceVertices  = mesh.faceVertices; 

  // Number of subcells
  Nsubcells = N*N;  
  // Number of nodes in this sub-triangulation
  Np = 0.5*(N+1)*(N+2);
  
  
// a very simple tesselation
if(settings.compareSetting("SUBCELL MINOR GRID","EQUISPACED")){  
  mesh.EquispacedNodesTri2D(N, vr, vs); 
  CellEquispacedEToVTri2D(N, mEToV);
}else if(settings.compareSetting("SUBCELL MINOR GRID","WARPBLEND")){
  mesh.NodesTri2D(N, vr, vs);   
  CellWarpBlendEToVTri2D(N, mEToV);
}

// Reference element subcell's center, face centers
cr.malloc(Nsubcells);
cs.malloc(Nsubcells);
mJ.malloc(Nsubcells);
// 
fr.malloc(Nsubcells*Nfaces);
fs.malloc(Nsubcells*Nfaces);

for(int s=0; s<Nsubcells; s++){
  dfloat tmpx = 0.0, tmpy = 0.0;
  // cell center 
  for(int v=0; v<Nverts; v++){
    const int vid = mEToV[s*Nverts+v];
    tmpx += vr[vid];   tmpy += vs[vid];
  }
  //
  cr[s] = tmpx/Nverts;  cs[s] = tmpy/Nverts;

  // face center
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

  // Cell Jacobian i.e. A_cell/A_element
  const int v1 = mEToV[s*Nverts+0]; 
  const int v2 = mEToV[s*Nverts+1]; 
  const int v3 = mEToV[s*Nverts+2]; 

  const dfloat xv1 = vr[v1], yv1 = vs[v1];
  const dfloat xv2 = vr[v2], yv2 = vs[v2];    
  const dfloat xv3 = vr[v3], yv3 = vs[v3];    
  mJ[s] = 0.5*( (xv2-xv1)*(yv3-yv1) - (xv3-xv1)*(yv2-yv1) )/2.0;
}


// Node of element(DG) to Cell at faces, required for mixed element lifting
// Not Sure yet AK....
mDGID.malloc(N*Nfaces); 
for(int n=0; n<Nfaces*N; n++){
 const int face        = n/N;  
 const int lid         = n%N;  
 if(lid<mesh.Nfp){
  const int dgid = face*N + lid; 
  mDGID[dgid]  = face*mesh.Nfp + lid; 
 }
}

o_mDGID = platform.malloc<dfloat>(mDGID); 

}


void stab_t::CellEquispacedEToVTri2D(const int _N, memory<int>& _EToV){
  const int _Nverts = 3;
  const int _Nelements = _N*_N;

  _EToV.malloc(_Nelements*_Nverts);

  int cnt=0;
  int sk=0;
  for (int j=0;j<_N;j++) {
    int shift = _N+1-j; //number of nodes in this row

    for (int i=0;i<_N-j;i++) { // upward triangle
      _EToV[cnt*_Nverts+0] = sk  ;
      _EToV[cnt*_Nverts+1] = sk+1;
      _EToV[cnt*_Nverts+2] = sk+shift;
      cnt++;
      if (i!=_N-j-1) { // downward triangle
        _EToV[cnt*_Nverts+0] = sk+shift+1; 
        _EToV[cnt*_Nverts+1] = sk+shift;
        _EToV[cnt*_Nverts+2] = sk+1;
        cnt++;
      }
      sk++;
    }
    sk++;
  }
}


void stab_t::CellWarpBlendEToVTri2D(const int _N, memory<int>& _EToV){
  const int _Nverts = 3;
  const int _Nelements = _N*_N;

  _EToV.malloc(_Nelements*_Nverts);

  int cnt=0;
  int sk=0;
  for (int j=0;j<_N;j++) {
    int shift = _N+1-j; //number of nodes in this row

    for (int i=0;i<_N-j;i++) { // upward triangle
      _EToV[cnt*_Nverts+0] = sk  ;
      _EToV[cnt*_Nverts+1] = sk+1;
      _EToV[cnt*_Nverts+2] = sk+shift;
      cnt++;
      if (i!=_N-j-1) { // downward triangle
        // _EToV[cnt*_Nverts+0] = sk+shift+1; 
        // _EToV[cnt*_Nverts+1] = sk+shift;
        // _EToV[cnt*_Nverts+2] = sk+1;
        _EToV[cnt*_Nverts+0] = sk+1;
        _EToV[cnt*_Nverts+1] = sk+shift+1; 
        _EToV[cnt*_Nverts+2] = sk+shift;

        cnt++;
      }
      sk++;
    }
    sk++;
  }

}







} //namespace libp
