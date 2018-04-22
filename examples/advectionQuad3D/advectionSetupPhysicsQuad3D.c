#include "advectionQuad3D.h"

void brownMinion(dfloat bmRho, dfloat bmDelta, dfloat sphereRadius,
		 dfloat x, dfloat y, dfloat z,
		 dfloat *u, dfloat *v, dfloat *w){

  dfloat Utangential = 0.25*(1+tanh(bmRho*(-z+0.5)))*(1+tanh(bmRho*(0.5+z)));

  dfloat uout, vout;

  if(x*x+y*y>1e-4) {
    uout =  -y*Utangential/(x*x+y*y);
    vout =   x*Utangential/(x*x+y*y);
  }
  else{
    uout = 0;
    vout = 0;
  }

  dfloat wout = bmDelta*sin(2*atan2(y,x))*(1-z*z);

  dfloat udotx = uout*x+vout*y+wout*z;
  *u = uout - udotx*x/(sphereRadius*sphereRadius);
  *v = vout - udotx*y/(sphereRadius*sphereRadius);
  *w = wout - udotx*z/(sphereRadius*sphereRadius);

}

solver_t *advectionSetupPhysicsQuad3D(mesh_t *mesh) {
  
  solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));
  solver->mesh = mesh;

  solver->Nfields = 10;

  // set temperature, gas constant, wave speeds
  solver->RT = 9.;
  solver->sqrtRT = sqrt(solver->RT);
  dfloat nu = 4e-4; //previous values: 8.9e-5; 9.8e-4; 6.e-3;
  
  // initial conditions
  dfloat sR = mesh->sphereRadius;

  solver->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*solver->Nfields,
				  sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      // Brown Minion shear layer roll up
      dfloat bmRho = 40; // was 40
      dfloat bmDelta  = 0.05;

      dfloat rho = 1;

      dfloat umod, vmod, wmod;

      brownMinion(bmRho, bmDelta, sR, x, y, z, &umod, &vmod, &wmod);

      dfloat delta = 1e-5;

      dfloat uP, uM, vP, vM, wP, wM;

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x+delta, y, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x-delta, y, z, &uM, &vM, &wM);

      dfloat dudx = (uP-uM)/(2*delta);
      dfloat dvdx = (vP-vM)/(2*delta);
      dfloat dwdx = (wP-wM)/(2*delta);

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y+delta, z, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y-delta, z, &uM, &vM, &wM);

      dfloat dudy = (uP-uM)/(2*delta);
      dfloat dvdy = (vP-vM)/(2*delta);
      dfloat dwdy = (wP-wM)/(2*delta);

      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z+delta, &uP, &vP, &wP);
      brownMinion(bmRho, bmDelta, mesh->sphereRadius, x, y, z-delta, &uM, &vM, &wM);

      dfloat dudz = (uP-uM)/(2*delta);
      dfloat dvdz = (vP-vM)/(2*delta);
      dfloat dwdz = (wP-wM)/(2*delta);

      dfloat divu = dudx + dvdy + dwdz;

#if 1
      dfloat sigma11 = nu*(dudx+dudx - (2*divu/3));
      dfloat sigma12 = nu*(dvdx+dudy);
      dfloat sigma13 = nu*(dwdx+dudz);
      dfloat sigma22 = nu*(dvdy+dvdy - (2*divu/3));
      dfloat sigma23 = nu*(dwdy+dvdz);
      dfloat sigma33 = nu*(dwdz+dwdz - (2*divu/3));
#else
      dfloat sigma11 = 0;
      dfloat sigma12 = 0;
      dfloat sigma13 = 0;
      dfloat sigma22 = 0;
      dfloat sigma23 = 0;
      dfloat sigma33 = 0;
#endif
      dfloat q1bar = rho;
      dfloat q2bar = rho*umod/mesh->sqrtRT;
      dfloat q3bar = rho*vmod/mesh->sqrtRT;
      dfloat q4bar = rho*wmod/mesh->sqrtRT;
      dfloat q5bar = (rho*umod*umod - sigma11)/(sqrt(2.)*mesh->RT);
      dfloat q6bar = (rho*vmod*vmod - sigma22)/(sqrt(2.)*mesh->RT);
      dfloat q7bar = (rho*wmod*wmod - sigma33)/(sqrt(2.)*mesh->RT);
      dfloat q8bar  = (rho*umod*vmod - sigma12)/mesh->RT;
      dfloat q9bar =  (rho*umod*wmod - sigma13)/mesh->RT;
      dfloat q10bar = (rho*vmod*wmod - sigma23)/mesh->RT;

      dfloat t = 0;

      int base = n + e*mesh->Np*solver->Nfields;

      /*      solver->q[base+0*mesh->Np] = q1bar; // uniform density, zero flow

	      solver->q[base+1*mesh->Np] = q2bar;
	      solver->q[base+2*mesh->Np] = q3bar;
	      solver->q[base+3*mesh->Np] = q4bar;

	      solver->q[base+4*mesh->Np] = q5bar;
	      solver->q[base+5*mesh->Np] = q6bar;
	      solver->q[base+6*mesh->Np] = q7bar;

	      solver->q[base+7*mesh->Np] = q8bar;
	      solver->q[base+8*mesh->Np] = q9bar;
	      solver->q[base+9*mesh->Np] = q10bar;*/

      solver->q[base+0*mesh->Np] = 1 + .1*exp(-20*((x-1)*(x-1)+y*y+z*z));
      solver->q[base+1*mesh->Np] = 1 ;
    }
  }

  // set BGK collision relaxation rate
  // nu = R*T*tau
  // 1/tau = RT/nu
  //  dfloat nu = 1.e-2/.5;
  //  dfloat nu = 1.e-3/.5;
  //  dfloat nu = 5.e-4;
  //    dfloat nu = 1.e-2; TW works for start up fence

  //set to 0 for advection
  solver->tauInv = 0;//mesh->RT/nu; // TW

  
  
  return solver;
}
