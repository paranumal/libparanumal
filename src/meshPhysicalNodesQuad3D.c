#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

//original mapping (r,s) -> (x,y)
void meshPhysicalNodesQuad3D(mesh_t *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      //      printf("xlin=%g, ylin=%g, zlin=%g\n", xlin, ylin, zlin);
      
      // project to sphere
      dfloat rlin = sqrt(xlin*xlin+ylin*ylin+zlin*zlin);
      mesh->x[cnt] = xlin/rlin; 
      mesh->y[cnt] = ylin/rlin; 
      mesh->z[cnt] = zlin/rlin; 
	
      ++cnt;
    }
  }
}

//original mapping (r,s) -> (x,y)
void meshSphericalNodesQuad3D(mesh_t *mesh){

  //constants used in conversions
  const dfloat R = 1;
  const dfloat a = 1./sqrt(3.);
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    int faceNumber = mesh->cubeFaceNumber[e];
    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of node on cube */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      dfloat xsph, ysph, zsph, cubRad;
      
      switch(faceNumber) {
      case 1: //positive constant x
	cubRad = sqrt(ylin*ylin + zlin*zlin + a*a);
	xsph = R/cubRad * a;
	ysph = R/cubRad * ylin;
	zsph = R/cubRad * zlin;
	break;
      case 2: //positive constant y
	cubRad = sqrt(xlin*xlin + zlin*zlin + a*a);
	xsph = R/cubRad * xlin;
	ysph = R/cubRad * a;
	zsph = R/cubRad * zlin;
	break;
      case 3: //negative constant x
	cubRad = sqrt(ylin*ylin + zlin*zlin + a*a);
	xsph =  -1 * R/cubRad * a;
	ysph =  R/cubRad * ylin;
	zsph =  R/cubRad * zlin;
	break;
      case 4: //negative constant y
	cubRad = sqrt(xlin*xlin + zlin*zlin + a*a);
	xsph = R/cubRad * xlin;
	ysph = -1*R/cubRad * a;
	zsph = R/cubRad * zlin;
	break;
      case 5: //positive constant z
	cubRad = sqrt(xlin*xlin + ylin*ylin + a*a);
	xsph = R/cubRad*xlin;
	ysph = R/cubRad*ylin;
	zsph = R/cubRad*a;
	break;
      case 6: //negative constant z
	cubRad = sqrt(xlin*xlin + ylin*ylin + a*a);
	xsph = R/cubRad * xlin;
	ysph = R/cubRad * ylin;
	zsph = -1*R/cubRad * a;
	break;
      }

      //Apply coordinate shift to vertex arrays
      /*if (rn == mesh->r[0] && sn == mesh->s[0]) {
	mesh->EX[id + 0] = xsph;
	mesh->EY[id + 0] = ysph;
	mesh->EZ[id + 0] = zsph;
      }
      else if (rn == mesh->r[0] && sn == mesh->s[mesh->Np - 1]) {
	mesh->EX[id + 1] = xsph;
	mesh->EY[id + 1] = ysph;
	mesh->EZ[id + 1] = zsph;
      }
      else if (rn == mesh->r[mesh->Np - 1] && sn == mesh->s[mesh->Np - 1]) {
	mesh->EX[id + 2] = xsph;
	mesh->EY[id + 2] = ysph;
	mesh->EZ[id + 2] = zsph;
      }
      else if (rn == mesh->r[0]&& sn == mesh->s[mesh->Np - 1]) {
	mesh->EX[id + 3] = xsph;
	mesh->EY[id + 3] = ysph;
	mesh->EZ[id + 3] = zsph;
	}*/
	  
      // project to sphere
      mesh->x[cnt] = xsph; 
      mesh->y[cnt] = ysph; 
      mesh->z[cnt] = zsph;
      
      ++cnt;
    }
  }
}

//untested effort to place spherical nodes on predetermined mesh
//Does not currently handle the poles (everything else should be present)
void meshSphericalNodesQuad3D_exp(mesh_t *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    //phi equatorial in the (x,y) plane, theta on the z axis
    //TODO: make sure acos doesn't break
    dfloat theta1 = acos(ze1/mesh->sphereRadius);
    dfloat theta2 = acos(ze2/mesh->sphereRadius);
    dfloat theta3 = acos(ze3/mesh->sphereRadius);
    dfloat theta4 = acos(ze4/mesh->sphereRadius);
    
    dfloat phi1 = atan2(ye1,xe1);
    dfloat phi2 = atan2(ye2,xe2);
    dfloat phi3 = atan2(ye3,xe3);
    dfloat phi4 = atan2(ye4,xe4);

    //correct azimuthal branch on nodes that straddle -M_PI

    //first go in a loop and make sure we take the short way around
    if (abs(phi2 - phi1) > M_PI) {
      if (phi2 > phi1) phi1 += M_PI;
      else phi2 += M_PI;
    }
    if (abs(phi3 - phi2) > M_PI) {
      if (phi3 > phi2) phi2 += M_PI;
      else phi3 += M_PI;
    }
    if (abs(phi4 - phi3) > M_PI) {
      if (phi4 > phi3) phi3 += M_PI;
      else phi4 += M_PI;
    }
    if (abs(phi1 - phi4) > M_PI) {
      if (phi1 > phi4) phi4 += M_PI;
      else phi1 += M_PI;
    }

    //some elements might have been stranded before, so check diagonals
    if (abs(phi1 - phi3) > M_PI) {
      if (phi1 > phi3) phi3 += M_PI;
      else phi1 += M_PI;
    }
    if (abs(phi2 - phi4) > M_PI) {
      if (phi2 > phi4) phi4 += M_PI;
      else phi2 += M_PI;
    }    

    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      dfloat thetalin = 
	+0.25*(1-rn)*(1-sn)*theta1
	+0.25*(1+rn)*(1-sn)*theta2
	+0.25*(1+rn)*(1+sn)*theta3
	+0.25*(1-rn)*(1+sn)*theta4;

      dfloat philin =
	+0.25*(1-rn)*(1-sn)*phi1
	+0.25*(1+rn)*(1-sn)*phi2
	+0.25*(1+rn)*(1+sn)*phi3
	+0.25*(1-rn)*(1+sn)*phi4;
      
      //      printf("xlin=%g, ylin=%g, zlin=%g\n", xlin, ylin, zlin);
      
      // project to sphere
      mesh->x[cnt] = mesh->sphereRadius * cos(philin) * sin(thetalin); 
      mesh->y[cnt] = mesh->sphereRadius * sin(philin) * sin(thetalin); 
      mesh->z[cnt] = mesh->sphereRadius * cos(thetalin); 
	
      ++cnt;
    }
  }
}
