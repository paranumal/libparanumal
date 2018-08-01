#include <mpi.h>

#include "fortranInterface.h"

mesh3D* meshParallelReaderHex3DExternal(MPI_Comm comm,
    hlong NHNnodes, hlong NHNhexes, hlong NHNboundaryFaces,
    hlong* EToV, hlong* BToV,
    dfloat* NHVX, dfloat* NHVY, dfloat* NHVZ) {

  int rank=0, size=1;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh3D *mesh = (mesh3D*) calloc(1, sizeof(mesh3D));
  MPI_Comm_dup(comm,&mesh->comm);

  mesh->dim = 3;
  mesh->rank = rank;
  mesh->size = size;
  mesh->Nverts = 8; // number of vertices per element
  mesh->Nfaces = 6;
  mesh->NfaceVertices = 4;

  // vertices on each face
  int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};

  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));

  memcpy(mesh->faceVertices, faceVertices[0],
         mesh->NfaceVertices*mesh->Nfaces*sizeof(int));

  /* read number of nodes in mesh */
  mesh->Nnodes = NHNnodes;

  dfloat *VX = NHVX;
  dfloat *VY = NHVY;
  dfloat *VZ = NHVZ;

  hlong Nelements = NHNhexes + NHNboundaryFaces;
  hlong Nhexes = 0, NboundaryFaces = NHNboundaryFaces;

  hlong NhexesLocal = NHNhexes;

  /* allocate space for Element node index data */
  mesh->EToV
    = (hlong*) calloc(NhexesLocal*mesh->Nverts, sizeof(hlong));

  mesh->elementInfo
    = (int*) calloc(NhexesLocal,sizeof(int));

  mesh->boundaryInfo = (hlong*) calloc(NboundaryFaces*
                                       (mesh->NfaceVertices+1), sizeof(hlong));

  for(hlong n=0; n<NboundaryFaces; ++n) {
    mesh->boundaryInfo[n*5+0] = BToV[n*(mesh->NfaceVertices+1)+0];
    mesh->boundaryInfo[n*5+1] = BToV[n*(mesh->NfaceVertices+1)+1]-1;
    mesh->boundaryInfo[n*5+2] = BToV[n*(mesh->NfaceVertices+1)+2]-1;
    mesh->boundaryInfo[n*5+3] = BToV[n*(mesh->NfaceVertices+1)+3]-1;
    mesh->boundaryInfo[n*5+4] = BToV[n*(mesh->NfaceVertices+1)+4]-1;
  }

  for(hlong n=0; n<NhexesLocal; ++n) {
    mesh->elementInfo[n] = 5;
    mesh->EToV[n*mesh->Nverts+0] = EToV[n*mesh->Nverts+0]-1;
    mesh->EToV[n*mesh->Nverts+1] = EToV[n*mesh->Nverts+1]-1;
    mesh->EToV[n*mesh->Nverts+2] = EToV[n*mesh->Nverts+2]-1;
    mesh->EToV[n*mesh->Nverts+3] = EToV[n*mesh->Nverts+3]-1;
    mesh->EToV[n*mesh->Nverts+4] = EToV[n*mesh->Nverts+4]-1;
    mesh->EToV[n*mesh->Nverts+5] = EToV[n*mesh->Nverts+5]-1;
    mesh->EToV[n*mesh->Nverts+6] = EToV[n*mesh->Nverts+6]-1;
    mesh->EToV[n*mesh->Nverts+7] = EToV[n*mesh->Nverts+7]-1;
  }

  /* record number of boundary faces found */
  mesh->NboundaryFaces = NboundaryFaces;

  /* record number of found hexes */
  mesh->Nelements = (hlong) NhexesLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts*mesh->Nelements,
                              sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts*mesh->Nelements,
                              sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts*mesh->Nelements,
                              sizeof(dfloat));
  for(hlong e=0; e<mesh->Nelements; ++e) {
    for(int n=0; n<mesh->Nverts; ++n) {
      hlong vid = mesh->EToV[e*mesh->Nverts+n];
      mesh->EX[e*mesh->Nverts+n] = VX[vid];
      mesh->EY[e*mesh->Nverts+n] = VY[vid];
      mesh->EZ[e*mesh->Nverts+n] = VZ[vid];
    }
  }

  return mesh;

}

void meshPhysicalNodesHex3DExternal(mesh_t *mesh, dfloat *x, dfloat *y, dfloat *z) {
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  hlong cnt = 0;
  for(int i=0; i<mesh->Nelements; i++) {
    for(int n=0; n<mesh->Np; n++) {
      mesh->x[cnt] = x[cnt];
      mesh->y[cnt] = y[cnt];
      mesh->z[cnt] = z[cnt];

      ++cnt;
    }
  }
}

void meshGeometricFactorsHex3DExternal(mesh3D *mesh,
  dfloat *NHxr, dfloat *NHxs, dfloat *NHxt,
  dfloat *NHyr, dfloat *NHys, dfloat *NHyt,
  dfloat *NHzr, dfloat *NHzs, dfloat *NHzt){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 12;
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  mesh->cubvgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->cubNp, sizeof(dfloat));

  /* number of second order geometric factors */
  mesh->Nggeo = 7;
  mesh->ggeo = (dfloat*) calloc(mesh->Nelements*mesh->Nggeo*mesh->Np, sizeof(dfloat));

  dfloat minJ = 1e9, maxJ = -1e9, maxSkew = 0;
  
  dbgfl;
  for(hlong e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    hlong id = e*mesh->Nverts;
    
    for(int k=0;k<mesh->Nq;++k){
      for(int j=0;j<mesh->Nq;++j){
        for(int i=0;i<mesh->Nq;++i){
          
	  int n = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq;
          hlong id = e*mesh->Np + n;
         
          dfloat xr = NHxr[id];
          dfloat xs = NHxs[id];
          dfloat xt = NHxt[id];
          dfloat yr = NHyr[id];
          dfloat ys = NHys[id];
          dfloat yt = NHyt[id];
          dfloat zr = NHzr[id];
          dfloat zs = NHzs[id];
          dfloat zt = NHzt[id];

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
          dfloat hr = sqrt(xr*xr+yr*yr+zr*zr);
          dfloat hs = sqrt(xs*xs+ys*ys+zs*zs);
          dfloat ht = sqrt(xt*xt+yt*yt+zt*zt);
          minJ = mymin(J, minJ);
          maxJ = mymax(J, maxJ);
          maxSkew = mymax(maxSkew, hr/hs);
          maxSkew = mymax(maxSkew, hr/ht);
          maxSkew = mymax(maxSkew, hs/hr);
          maxSkew = mymax(maxSkew, hs/ht);
          maxSkew = mymax(maxSkew, ht/hr);
          maxSkew = mymax(maxSkew, ht/hs);
	  //  dbgfl;
          
          if(J<1e-12) printf("J = %g !!!!!!!!!!!!!\n", J);
          
          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
          
          dfloat JW = J*mesh->gllw[i]*mesh->gllw[j]*mesh->gllw[k];
          // dbgfl;
          
          /* store geometric factors */
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RXID] = rx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RYID] = ry;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*RZID] = rz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SXID] = sx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SYID] = sy;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*SZID] = sz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TXID] = tx;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TYID] = ty;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*TZID] = tz;
          
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JID]  = J;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*JWID] = JW;
          mesh->vgeo[mesh->Nvgeo*mesh->Np*e + n + mesh->Np*IJWID] = 1./JW;

          /* store second order geometric factors */
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
          mesh->ggeo[mesh->Nggeo*mesh->Np*e + n + mesh->Np*GWJID] = JW;
        }
      }
    }

    //TODO
    //geometric data for quadrature
//  for(int k=0;k<mesh->cubNq;++k){
//    for(int j=0;j<mesh->cubNq;++j){
//      for(int i=0;i<mesh->cubNq;++i){
//        
//        dfloat rn = mesh->cubr[i];
//        dfloat sn = mesh->cubr[j];
//        dfloat tn = mesh->cubr[k];

//        /* Jacobian matrix */
//        dfloat xr = 0.125*( (1-tn)*(1-sn)*(xe[1]-xe[0]) + (1-tn)*(1+sn)*(xe[2]-xe[3]) + (1+tn)*(1-sn)*(xe[5]-xe[4]) + (1+tn)*(1+sn)*(xe[6]-xe[7]) );
//        dfloat xs = 0.125*( (1-tn)*(1-rn)*(xe[3]-xe[0]) + (1-tn)*(1+rn)*(xe[2]-xe[1]) + (1+tn)*(1-rn)*(xe[7]-xe[4]) + (1+tn)*(1+rn)*(xe[6]-xe[5]) );
//        dfloat xt = 0.125*( (1-rn)*(1-sn)*(xe[4]-xe[0]) + (1+rn)*(1-sn)*(xe[5]-xe[1]) + (1+rn)*(1+sn)*(xe[6]-xe[2]) + (1-rn)*(1+sn)*(xe[7]-xe[3]) );
//        
//        dfloat yr = 0.125*( (1-tn)*(1-sn)*(ye[1]-ye[0]) + (1-tn)*(1+sn)*(ye[2]-ye[3]) + (1+tn)*(1-sn)*(ye[5]-ye[4]) + (1+tn)*(1+sn)*(ye[6]-ye[7]) );
//        dfloat ys = 0.125*( (1-tn)*(1-rn)*(ye[3]-ye[0]) + (1-tn)*(1+rn)*(ye[2]-ye[1]) + (1+tn)*(1-rn)*(ye[7]-ye[4]) + (1+tn)*(1+rn)*(ye[6]-ye[5]) );
//        dfloat yt = 0.125*( (1-rn)*(1-sn)*(ye[4]-ye[0]) + (1+rn)*(1-sn)*(ye[5]-ye[1]) + (1+rn)*(1+sn)*(ye[6]-ye[2]) + (1-rn)*(1+sn)*(ye[7]-ye[3]) );
//        
//        dfloat zr = 0.125*( (1-tn)*(1-sn)*(ze[1]-ze[0]) + (1-tn)*(1+sn)*(ze[2]-ze[3]) + (1+tn)*(1-sn)*(ze[5]-ze[4]) + (1+tn)*(1+sn)*(ze[6]-ze[7]) );
//        dfloat zs = 0.125*( (1-tn)*(1-rn)*(ze[3]-ze[0]) + (1-tn)*(1+rn)*(ze[2]-ze[1]) + (1+tn)*(1-rn)*(ze[7]-ze[4]) + (1+tn)*(1+rn)*(ze[6]-ze[5]) );
//        dfloat zt = 0.125*( (1-rn)*(1-sn)*(ze[4]-ze[0]) + (1+rn)*(1-sn)*(ze[5]-ze[1]) + (1+rn)*(1+sn)*(ze[6]-ze[2]) + (1-rn)*(1+sn)*(ze[7]-ze[3]) );
//        
//        /* compute geometric factors for affine coordinate transform*/
//        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
//        
//        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
//        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
//        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
//        
//        dfloat JW = J*mesh->cubw[i]*mesh->cubw[j]*mesh->cubw[k];
//        
//        /* store geometric factors */
//        dlong base = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq + k*mesh->cubNq*mesh->cubNq;
//        mesh->cubvgeo[base + mesh->cubNp*RXID] = rx;
//        mesh->cubvgeo[base + mesh->cubNp*RYID] = ry;
//        mesh->cubvgeo[base + mesh->cubNp*RZID] = rz;
//        
//        mesh->cubvgeo[base + mesh->cubNp*SXID] = sx;
//        mesh->cubvgeo[base + mesh->cubNp*SYID] = sy;
//        mesh->cubvgeo[base + mesh->cubNp*SZID] = sz;
//        
//        mesh->cubvgeo[base + mesh->cubNp*TXID] = tx;
//        mesh->cubvgeo[base + mesh->cubNp*TYID] = ty;
//        mesh->cubvgeo[base + mesh->cubNp*TZID] = tz;
//        
//        mesh->cubvgeo[base + mesh->cubNp*JID]  = J;
//        mesh->cubvgeo[base + mesh->cubNp*JWID] = JW;
//        mesh->cubvgeo[base + mesh->cubNp*IJWID] = 1./JW;
//      }
//    }
//  }
  }

  printf("J in range [%g,%g] and max Skew = %g\n", minJ, maxJ, maxSkew);
}


/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsHex3DExternal(mesh3D *mesh,
  dfloat *NHxr, dfloat *NHxs, dfloat *NHxt,
  dfloat *NHyr, dfloat *NHys, dfloat *NHyt,
  dfloat *NHzr, dfloat *NHzs, dfloat *NHzt){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 14;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
                                sizeof(dfloat));

  mesh->cubsgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
                                mesh->Nsgeo*mesh->cubNfp*mesh->Nfaces, 
                                sizeof(dfloat));

  for(hlong e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    hlong id = e*mesh->Nverts;

    for(int f=0;f<mesh->Nfaces;++f){ // for each face
      
      for(int i=0;i<mesh->Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = mesh->faceNodes[f*mesh->Nfp+i];
	hlong id = e*mesh->Np+n;

        dfloat xr = NHxr[id];
        dfloat xs = NHxs[id];
        dfloat xt = NHxt[id];
        dfloat yr = NHyr[id];
        dfloat ys = NHys[id];
        dfloat yt = NHyt[id];
        dfloat zr = NHzr[id];
        dfloat zs = NHzs[id];
        dfloat zt = NHzt[id];

        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
        
        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
        
        /* face f normal and length */
        dfloat nx, ny, nz, d;
        switch(f){
        case 0: nx = -tx; ny = -ty; nz = -tz; break;
        case 1: nx = -sx; ny = -sy; nz = -sz; break;
        case 2: nx = +rx; ny = +ry; nz = +rz; break;
        case 3: nx = +sx; ny = +sy; nz = +sz; break;
        case 4: nx = -rx; ny = -ry; nz = -rz; break;
        case 5: nx = +tx; ny = +ty; nz = +tz; break;
        }

        dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
        nx /= sJ; ny /= sJ; nz /= sJ;
        sJ *= J;
        
        /* output index */
        hlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->Nfp*e + mesh->Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->sgeo[base+NXID] = nx;
        mesh->sgeo[base+NYID] = ny;
        mesh->sgeo[base+NZID] = nz;
        mesh->sgeo[base+SJID] = sJ;
        mesh->sgeo[base+IJID] = 1./J;

        mesh->sgeo[base+WIJID] = 1./(J*mesh->gllw[0]);
        mesh->sgeo[base+WSJID] = sJ*mesh->gllw[i%mesh->Nq]*mesh->gllw[i/mesh->Nq];

	computeFrame(nx, ny, nz,
		     mesh->sgeo[base+STXID], mesh->sgeo[base+STYID], mesh->sgeo[base+STZID],
		     mesh->sgeo[base+SBXID], mesh->sgeo[base+SBYID], mesh->sgeo[base+SBZID]);
      }

      //TODO
      //geometric data for quadrature
//    for(int i=0;i<mesh->cubNfp;++i){  // for each quadrature node on face

//      dfloat rn, sn, tn;
//      switch(f){
//      case 0: rn = mesh->cubr[i%mesh->cubNq]; sn = mesh->cubr[i/mesh->cubNq]; tn = -1.0;                      break;
//      case 1: rn = mesh->cubr[i%mesh->cubNq]; sn = -1.0;                      tn = mesh->cubr[i/mesh->cubNq]; break;
//      case 2: rn = 1.0;                       sn = mesh->cubr[i%mesh->cubNq]; tn = mesh->cubr[i/mesh->cubNq]; break;
//      case 3: rn = mesh->cubr[i%mesh->cubNq]; sn = 1.0;                       tn = mesh->cubr[i/mesh->cubNq]; break;
//      case 4: rn = -1.0;                      sn = mesh->cubr[i%mesh->cubNq]; tn = mesh->cubr[i/mesh->cubNq]; break;
//      case 5: rn = mesh->cubr[i%mesh->cubNq]; sn = mesh->cubr[i/mesh->cubNq]; tn = 1.0;                       break;
//      }

//      /* Jacobian matrix */
//      dfloat xr = 0.125*( (1-tn)*(1-sn)*(xe[1]-xe[0]) + (1-tn)*(1+sn)*(xe[2]-xe[3]) + (1+tn)*(1-sn)*(xe[5]-xe[4]) + (1+tn)*(1+sn)*(xe[6]-xe[7]) );
//      dfloat xs = 0.125*( (1-tn)*(1-rn)*(xe[3]-xe[0]) + (1-tn)*(1+rn)*(xe[2]-xe[1]) + (1+tn)*(1-rn)*(xe[7]-xe[4]) + (1+tn)*(1+rn)*(xe[6]-xe[5]) );
//      dfloat xt = 0.125*( (1-rn)*(1-sn)*(xe[4]-xe[0]) + (1+rn)*(1-sn)*(xe[5]-xe[1]) + (1+rn)*(1+sn)*(xe[6]-xe[2]) + (1-rn)*(1+sn)*(xe[7]-xe[3]) );
//      
//      dfloat yr = 0.125*( (1-tn)*(1-sn)*(ye[1]-ye[0]) + (1-tn)*(1+sn)*(ye[2]-ye[3]) + (1+tn)*(1-sn)*(ye[5]-ye[4]) + (1+tn)*(1+sn)*(ye[6]-ye[7]) );
//      dfloat ys = 0.125*( (1-tn)*(1-rn)*(ye[3]-ye[0]) + (1-tn)*(1+rn)*(ye[2]-ye[1]) + (1+tn)*(1-rn)*(ye[7]-ye[4]) + (1+tn)*(1+rn)*(ye[6]-ye[5]) );
//      dfloat yt = 0.125*( (1-rn)*(1-sn)*(ye[4]-ye[0]) + (1+rn)*(1-sn)*(ye[5]-ye[1]) + (1+rn)*(1+sn)*(ye[6]-ye[2]) + (1-rn)*(1+sn)*(ye[7]-ye[3]) );
//      
//      dfloat zr = 0.125*( (1-tn)*(1-sn)*(ze[1]-ze[0]) + (1-tn)*(1+sn)*(ze[2]-ze[3]) + (1+tn)*(1-sn)*(ze[5]-ze[4]) + (1+tn)*(1+sn)*(ze[6]-ze[7]) );
//      dfloat zs = 0.125*( (1-tn)*(1-rn)*(ze[3]-ze[0]) + (1-tn)*(1+rn)*(ze[2]-ze[1]) + (1+tn)*(1-rn)*(ze[7]-ze[4]) + (1+tn)*(1+rn)*(ze[6]-ze[5]) );
//      dfloat zt = 0.125*( (1-rn)*(1-sn)*(ze[4]-ze[0]) + (1+rn)*(1-sn)*(ze[5]-ze[1]) + (1+rn)*(1+sn)*(ze[6]-ze[2]) + (1-rn)*(1+sn)*(ze[7]-ze[3]) );

//      /* determinant of Jacobian matrix */
//      dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
//      
//      dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
//      dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
//      dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
//      
//      /* face f normal and length */
//      dfloat nx, ny, nz, d;
//      switch(f){
//      case 0: nx = -tx; ny = -ty; nz = -tz; break;
//      case 1: nx = -sx; ny = -sy; nz = -sz; break;
//      case 2: nx = +rx; ny = +ry; nz = +rz; break;
//      case 3: nx = +sx; ny = +sy; nz = +sz; break;
//      case 4: nx = -rx; ny = -ry; nz = -rz; break;
//      case 5: nx = +tx; ny = +ty; nz = +tz; break;
//      }

//      dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
//      nx /= sJ; ny /= sJ; nz /= sJ;
//      sJ *= J;
//      

//      /* output index */
//      hlong base = mesh->Nsgeo*(mesh->Nfaces*mesh->cubNfp*e + mesh->cubNfp*f + i);

//      /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
//      mesh->cubsgeo[base+NXID] = nx;
//      mesh->cubsgeo[base+NYID] = ny;
//      mesh->cubsgeo[base+NZID] = nz;
//      mesh->cubsgeo[base+SJID] = sJ;
//      mesh->cubsgeo[base+IJID] = 1./J;

//      mesh->cubsgeo[base+WIJID] = 1./(J*mesh->cubw[0]);
//      mesh->cubsgeo[base+WSJID] = sJ*mesh->cubw[i%mesh->cubNq]*mesh->cubw[i/mesh->cubNq];

//      computeFrame(nx, ny, nz,
//      	     mesh->cubsgeo[base+STXID], mesh->cubsgeo[base+STYID], mesh->cubsgeo[base+STZID],
//      	     mesh->cubsgeo[base+SBXID], mesh->cubsgeo[base+SBYID], mesh->cubsgeo[base+SBZID]);
//      
//      
//    }
    }
  }

  for(hlong e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      hlong baseM = e*mesh->Nfp*mesh->Nfaces + n;
      hlong baseP = mesh->mapP[baseM];
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }
}
mesh_t *meshSetupHex3DExternal(mesh_t *mesh, int N,
  dfloat *x   , dfloat *y   , dfloat *z   ,
  dfloat *NHxr, dfloat *NHxs, dfloat *NHxt,
  dfloat *NHyr, dfloat *NHys, dfloat *NHyt,
  dfloat *NHzr, dfloat *NHzs, dfloat *NHzt){

  // @NEK-HOLMES: Commented out because we do this through Fortran-API
  // read chunk of elements
//mesh3D *mesh = meshParallelReaderHex3D(filename);

  // partition elements using Morton ordering & parallel sort
//  meshGeometricPartition3D(mesh); 
  
  // connect elements using parallel sort
  meshParallelConnect(mesh);
  
  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesHex3DExternal(mesh,x,y,z);

  // compute geometric factors

  meshGeometricFactorsHex3DExternal(mesh,
    NHxr, NHxs, NHxt, NHyr, NHys, NHyt,
    NHzr, NHzs, NHzt);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  meshConnectFaceNodes3D(mesh);

  meshSurfaceGeometricFactorsHex3DExternal(mesh,
    NHxr, NHxs, NHxt, NHyr, NHys, NHyt,
    NHzr, NHzs, NHzt);

  // global nodes
  meshParallelConnectNodes(mesh);

  // initialize LSERK4 time stepping coefficients
  int Nrk = 5;

  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  dfloat rkc[6] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0,
		   1.}; 

  mesh->Nrk = Nrk;
  memcpy(mesh->rka, rka, Nrk*sizeof(dfloat));
  memcpy(mesh->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(mesh->rkc, rkc, (Nrk+1)*sizeof(dfloat));
    
  return mesh;
}
