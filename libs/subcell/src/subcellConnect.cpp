#include "subcell.hpp"

// structure used to encode vertices that make
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  dlong element;
  int face;

  dlong elementNeighbor; // neighbor element
  int faceNeighbor;    // neighbor face

  int NfaceVertices;

  hlong v[4];

}face_t;


static void findBestMatch(dfloat x1, dfloat y1, int Nfaces,
			  int N, int *eList, dfloat *x2, dfloat *y2, int *nE, int *nF){

  int eIndex = 0;
  int fIndex = 0;
  dfloat mindist2 = 1e12;   

  for(int e=0;e<N;++e){
    /* next element */
    const int e2 = eList[e];
    /* check faces */
    for(int f2=0; f2<Nfaces; f2++){  
      /* distance between target and next node */
      const dfloat dist2 = pow(x1-x2[e2*Nfaces + f2],2) + 
	pow(y1-y2[e2*Nfaces + f2],2);
      /* if next node is closer to target update match */
      if(dist2<mindist2){
	mindist2 = dist2;
	eIndex = e2; fIndex = f2;
      }

    }
  }
   
  if(mindist2>1e-3) {
    stringstream ss;
    ss << "Bad match: x,y = " << x1 << ", " << y1 << "\n";
    LIBP_ABORT(ss.str())
      }
  
  *nE = eIndex; 
  *nF = fIndex; 
}

// comparison function that orders vertices
// based on their combined vertex indices
static int compareVertices(const void *a,
			   const void *b){

  face_t *fa = (face_t*) a;
  face_t *fb = (face_t*) b;

  for(int n=0;n<fa->NfaceVertices;++n){
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;

}
/* comparison function that orders element/face
   based on their indexes */
static int compareFaces(const void *a,
			const void *b){

  face_t *fa = (face_t*) a;
  face_t *fb = (face_t*) b;

  if(fa->element < fb->element) return -1;
  if(fa->element > fb->element) return +1;

  if(fa->face < fb->face) return -1;
  if(fa->face > fb->face) return +1;

  return 0;

}

/* routine to find EToE (Element To Element)
   and EToF (Element To Local Face) connectivity arrays */
void subcell_t::LocalConnect(){
  /* build list of faces */
  face_t *faces =
    (face_t*) calloc(Nsubcells*Nfaces, sizeof(face_t));

  dlong cnt = 0;
  for(dlong e=0;e<Nsubcells;++e){
    for(int f=0;f<Nfaces;++f){

      for(int n=0;n<NfaceVertices;++n){
        dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
        faces[cnt].v[n] = mEToV[vid];
      }

      mysort(faces[cnt].v, NfaceVertices, "descending");

      faces[cnt].NfaceVertices = NfaceVertices;

      faces[cnt].element = e;
      faces[cnt].face = f;

      faces[cnt].elementNeighbor= -1;
      faces[cnt].faceNeighbor = -1;

      ++cnt;
    }
  }

  /* sort faces by their vertex number pairs */
  qsort(faces, Nsubcells*Nfaces, sizeof(face_t), compareVertices);

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for(cnt=0;cnt<Nsubcells*Nfaces-1;++cnt){
    if(!compareVertices(faces+cnt, faces+cnt+1)){
      // match
      faces[cnt].elementNeighbor = faces[cnt+1].element;
      faces[cnt].faceNeighbor = faces[cnt+1].face;

      faces[cnt+1].elementNeighbor = faces[cnt].element;
      faces[cnt+1].faceNeighbor = faces[cnt].face;
    }
  }

  /* resort faces back to the original element/face ordering */
  qsort(faces, Nsubcells*Nfaces, sizeof(face_t), compareFaces);

  /* extract the element to element and element to face connectivity */
  mEToE = (int*) calloc(Nsubcells*Nfaces, sizeof(int));
  mEToF = (int*) calloc(Nsubcells*Nfaces, sizeof(int));

  cnt = 0;
  for(dlong e=0;e<Nsubcells;++e){
    for(int f=0;f<Nfaces;++f){
      mEToE[cnt] = faces[cnt].elementNeighbor;
      mEToF[cnt] = faces[cnt].faceNeighbor;
      ++cnt;
    }
  }

  // external connection 
  dfloat deps = 1.;
  while((1.+deps)>1.) deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;
  
  int fcnt[3]; 
  for (int i=0;i<3;i++) fcnt[i]=0;

  // local cell to subcell trace info   
  mFToE = (int*) calloc(N*Nfaces, sizeof(int));
  mFToF = (int*) calloc(N*Nfaces, sizeof(int));

  for(dlong e=0;e<Nsubcells;++e){
    for(int f=0;f<Nfaces;++f){
      if(mEToE[e*Nfaces + f]<0){ // local boundary element

        dfloat rf = 0, sf = 0; 
        for(int n=0;n<NfaceVertices;++n){
          dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
          rf += vr[mEToV[vid]]/NfaceVertices; 
          sf += vs[mEToV[vid]]/NfaceVertices;
        }

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

/* routine to find EToE (Element To Element)
   and EToF (Element To Local Face) connectivity arrays */
void subcell_t::GlobalConnect(){
  const dlong Nelements = mesh.Nelements; 
  emapP = (dlong*) calloc(Nsubcells*Nfaces*Nelements, sizeof(dlong));
  fmapP = (dlong*) calloc(Nsubcells*Nfaces*Nelements, sizeof(dlong));
  // first connect internal elements 
  for(dlong e=0; e<Nelements; e++){
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

  //check if we're using a periodic box mesh
  int periodicFlag = 0;
  if (mesh.settings.compareSetting("MESH FILE","BOX") &&
      mesh.settings.compareSetting("BOX BOUNDARY FLAG","-1"))
    periodicFlag = 1;

  //box dimensions
  dfloat DIMX, DIMY;
  mesh.settings.getSetting("BOX DIMX", DIMX);
  mesh.settings.getSetting("BOX DIMY", DIMY);

  //box is centered at the origin
  DIMX /= 2.0;
  DIMY /= 2.0;

  // Compute Face Centers and Connect
  dfloat *xf = (dfloat *) malloc(mesh.Nelements*Nsubcells*Nfaces*sizeof(dfloat));
  dfloat *yf = (dfloat *) malloc(mesh.Nelements*Nsubcells*Nfaces*sizeof(dfloat));

  dlong cnt = 0; 
  for(dlong e=0; e<mesh.Nelements; e++){
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

  for(dlong e=0; e<Nelements; e++){
    for(int f=0;f<Nfaces;++f){
      dlong eP = mesh.EToE[e*Nfaces + f]; 
      int fP   = mesh.EToF[e*Nfaces + f]; 

      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e; fP = f;
      }

      dfloat offsetX = 0.0;
      dfloat offsetY = 0.0;

      if (periodicFlag) {
        //if the mesh is periodic, this is more complicated.
        // check if this face is on a boundary face
        bool top=true, bottom=true, left=true, right=true;
        for(int n=0;n<NfaceVertices;++n){
          dlong vid = e*mesh.Nverts + mesh.faceVertices[f*NfaceVertices+n];
          if (fabs(mesh.EX[vid]-DIMX)>1e-4) right = false;
          if (fabs(mesh.EX[vid]+DIMX)>1e-4) left = false;
          if (fabs(mesh.EY[vid]-DIMY)>1e-4) top = false;
          if (fabs(mesh.EY[vid]+DIMY)>1e-4) bottom = false;
        }
        if (right)  offsetX = -2.0*DIMX;
        if (left)   offsetX =  2.0*DIMX;
        if (top)    offsetY = -2.0*DIMY;
        if (bottom) offsetY =  2.0*DIMY;
      }

      // for each subcell at this face find neighboor element
      for(int n=0; n<N; n++){
	// local element and face info 
	const int sem =  mFToE[f*N+ n];
	const int sfm =  mFToF[f*N+ n];
	//
	const dlong idM = e*Nsubcells*Nfaces + sem*Nfaces + sfm; 
        dfloat xM = xf[idM] + offsetX;
        dfloat yM = yf[idM] + offsetY;
        
        int idE, idF; // mEToE[sem*Nfaces + sfm]; 
        findBestMatch(xM, yM, Nfaces, N, mFToE+fP*N,
		      xf+eP*Nsubcells*Nfaces,
		      yf+eP*Nsubcells*Nfaces, &idE, &idF);

	const dlong eshift = e*Nsubcells*Nfaces + sem*Nfaces + sfm;  
	emapP[eshift]     = idE + eP*Nsubcells; 
	fmapP[eshift]     = idF; //

      }
    }
  }

  free(xf); 
  free(yf);
  
#if 0
  // first connect internal elements 
  for(dlong e=0; e<Nelements; e++){
    for(int s =0; s<Nsubcells; s++){
      for(int f=0; f<Nfaces; f++){
	const dlong sfm = e*Nsubcells*Nfaces + s*Nfaces + f; 
	const dlong sfp = emapP[sfm]*Nfaces + fmapP[sfm]; 
	if( fabs( x[sfm] - x[sfp] ) > 1e-12){
          stringstream ss;
          ss << "bad index= " << e << ", " << s << ", "<<f<<"\n";
          LIBP_ABORT(ss.str())
	    }
      }
    }
  }
#endif

}
