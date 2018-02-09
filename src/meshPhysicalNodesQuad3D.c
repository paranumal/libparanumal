#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

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

void meshSphericalNodesQuad3D(mesh_t *mesh){
  
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
      if (phi2 > phi1) phi1 += 2*M_PI;
      else phi2 += 2*M_PI;
    }
    if (abs(phi3 - phi2) > M_PI) {
      if (phi3 > phi2) phi2 += 2*M_PI;
      else phi3 += 2*M_PI;
    }
    if (abs(phi4 - phi3) > M_PI) {
      if (phi4 > phi3) phi3 += 2*M_PI;
      else phi4 += 2*M_PI;
    }
    if (abs(phi1 - phi4) > M_PI) {
      if (phi1 > phi4) phi4 += 2*M_PI;
      else phi1 += 2*M_PI;
    }

    //some elements might have been stranded before, so check diagonals
    if (abs(phi1 - phi2) > M_PI) {
      if (phi1 > phi2) phi2 += 2*M_PI;
      else phi1 += 2*M_PI;
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
