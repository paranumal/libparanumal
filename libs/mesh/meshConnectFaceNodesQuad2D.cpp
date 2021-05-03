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

#include "mesh.hpp"
#include "mesh/mesh2D.hpp"
#include "mesh/mesh3D.hpp"

void meshQuad3D::ConnectFaceNodes(){
  mesh_t *mesh_p = (mesh_t*) this;
  meshQuad2D* quadmesh = (meshQuad2D*) mesh_p;
  quadmesh->meshQuad2D::ConnectFaceNodes();
}

// serial face-node to face-node connection
void meshQuad2D::ConnectFaceNodes(){

  const dfloat NODETOL = 1.0e-5;

  dfloat V0[2] = {-1.0, 1.0};
  dfloat V1[2] = { 1.0,-1.0};

  dfloat EX0[Nverts];
  dfloat EX1[Nverts];

  dfloat *x0 = (dfloat*) malloc(Nfp*sizeof(dfloat));
  dfloat *x1 = (dfloat*) malloc(Nfp*sizeof(dfloat));

  /* Build the permutation array R */
  int *R = (int*) malloc(Nfaces*Nfaces*Nfp*sizeof(int));

  for (int fM=0;fM<Nfaces;fM++) {

    for (int v=0;v<Nverts;v++) EX0[v] = 0.0;

    //setup top element with face fM on the bottom
    for (int v=0;v<NfaceVertices;v++) {
      int fv = faceVertices[fM*NfaceVertices + v];
      EX0[fv] = V0[v];
    }

    for(int n=0;n<Nfp;++n){ /* for each face node */
      const int fn = faceNodes[fM*Nfp+n];

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = r[fn];
      dfloat sn = s[fn];

      /* physical coordinate of interpolation node */
      x0[n] =  0.25*(1-rn)*(1-sn)*EX0[0]
              +0.25*(1+rn)*(1-sn)*EX0[1]
              +0.25*(1+rn)*(1+sn)*EX0[2]
              +0.25*(1-rn)*(1+sn)*EX0[3];
    }

    for (int fP=0;fP<Nfaces;fP++) { /*For each neighbor face */
      // Zero vertices
      for (int v=0;v<Nverts;v++) EX1[v] = 0.0;

      //setup bottom element with face fP on the top
      for (int v=0;v<NfaceVertices;v++) {
        int fv = faceVertices[fP*NfaceVertices + v];
        EX1[fv] = V1[v];
      }

      for(int n=0;n<Nfp;++n){ /* for each node */
        const int fn = faceNodes[fP*Nfp+n];

        /* (r,s) coordinates of interpolation nodes*/
        dfloat rn = r[fn];
        dfloat sn = s[fn];

        /* physical coordinate of interpolation node */
        x1[n] =  0.25*(1-rn)*(1-sn)*EX1[0]
                +0.25*(1+rn)*(1-sn)*EX1[1]
                +0.25*(1+rn)*(1+sn)*EX1[2]
                +0.25*(1-rn)*(1+sn)*EX1[3];
      }

      /* for each node on this face find the neighbor node */
      for(int n=0;n<Nfp;++n){
        const dfloat xM = x0[n];

        int m=0;
        for(;m<Nfp;++m){ /* for each neighbor node */
          const dfloat xP = x1[m];

          /* distance between target and neighbor node */
          const dfloat dist = fabs(xM-xP);

          /* if neighbor node is close to target, match */
          if(dist<NODETOL){
            R[fM*Nfaces*Nfp + fP*Nfp + n] = m;
            break;
          }
        }

        /*Check*/
        const dfloat xP = x1[m];

        /* distance between target and neighbor node */
        const dfloat dist = fabs(xM-xP);
        if(dist>NODETOL){
          //This shouldn't happen
          stringstream ss;
          ss << "Unable to match face node, face: " << fM
             << ", matching face: " << fP
             << ", node: " << n
             << ". Is the reference node set not symmetric?";
          LIBP_ABORT(ss.str())
        }
      }
    }
  }

  free(x0);
  free(x1);

  /* volume indices of the interior and exterior face nodes for each element */
  vmapM = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  vmapP = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  mapP  = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));

  /* assume elements already connected */
  for(dlong eM=0;eM<Nelements;++eM){
    for(int fM=0;fM<Nfaces;++fM){
      dlong eP = EToE[eM*Nfaces+fM];
      int fP = EToF[eM*Nfaces+fM];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        for(int nM=0;nM<Nfp;++nM){
          const int idM = faceNodes[fM*Nfp+nM];
          const dlong id = eM*Nfaces*Nfp + fM*Nfp + nM;
          vmapM[id] = idM + eM*Np;
          vmapP[id] = idM + eM*Np;
          mapP[id] = eM*Nfaces*Nfp + fM*Nfp + nM;
        }
      } else {
        /* for each node on this face use the permuation array
           to select the neighbor node */
        for(int nM=0;nM<Nfp;++nM){
          const int nP  = R[fM*Nfaces*Nfp + fP*Nfp + nM];
          const int idM = faceNodes[fM*Nfp+nM];
          const int idP = faceNodes[fP*Nfp+nP];

          const dlong id = eM*Nfaces*Nfp + fM*Nfp + nM;
          vmapM[id] = idM + eM*Np;
          vmapP[id] = idP + eP*Np;
          mapP[id] = eP*Nfaces*Nfp + fP*Nfp + nP;

        #if 0
          /*Sanity check*/
          dfloat xnM = x[idM + eM*Np];
          dfloat ynM = y[idM + eM*Np];
          dfloat xnP = x[idP + eP*Np];
          dfloat ynP = y[idP + eP*Np];
          const dfloat dist = pow(xnM-xnP,2) + pow(ynM-ynP,2);
          if (dist>NODETOL)
            printf("Mismatch?: Element %d, face %d, node %d, xM = (%f, %f), xP = (%f, %f)\n",
                    eM, fM, nM, xnM, ynM, xnP, ynP);
        #endif
        }
      }
    }
  }

//      printf("connecting (%d,%d) to (%d,%d) [ vmapM %d to vmapP %d ]\n",
//             e,f,eP,fP, vmapM[id], vmapP[id]);
  free(R);
}

