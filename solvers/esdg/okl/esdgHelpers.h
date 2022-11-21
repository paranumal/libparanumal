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
#define p_gamma (1.4)
#define p_Pr 0.63

#define mymin(a,b) ( ((a)<(b)) ? (a):(b) )
#define mymax(a,b) ( ((a)>(b)) ? (a):(b) )

#ifndef GUARD_LEVEL
#define GUARD_LEVEL 2
#endif

// GUARD_LEVEL bit flags:
// 1 => guard sqrt  (sqrt(|r|))
// 2 => guard abs   (|r|) 
// 4 => guard log   (log(|r|))
// 8 => guard a/b   (a/(b+GUARD_DENOMINATOR))
//16 => guard pow(a,b)   pow(fabs(a),b)
dfloat guardSqrt(dfloat r){
#if GUARD_LEVEL&1
  return sqrt(fabs(r));
#else
  return sqrt(r);
#endif
}

dfloat guardPositive(dfloat r){
#if GUARD_LEVEL&2
  return fabs(r);
#else
  return r;
#endif
}

dfloat guardLog(dfloat r){
#if GUARD_LEVEL&4
  return log(fabs(r));
#else
  return log(r);
#endif
}

// assumes denominator (b) should be positive, bounded away from zero
#define GUARD_DENOMINATOR 1e-5
dfloat guardDivide(dfloat a, dfloat b){
#if GUARD_LEVEL&8
  return a/(fabs(b)+GUARD_DENOMINATOR);
#else
  return a/b;
#endif
}

// use when raising number to power that should be positive
dfloat guardPow(dfloat a, dfloat b){
#if GUARD_LEVEL&16
  return pow(fabs(a),b);
#else
  return pow(a,b);
#endif
}


#define PSCALE  (0.5)
//#define PSCALE  (1.0)
//#define PSCALE  (0.0)
// 0 seems to be broken
#define USE_LLF_PENALTY 1
#define USE_HLL_PENALTY 0
// Prandtl number (O2)
#define PRINT_WARNINGS 1

// flux forms from: https://arxiv.org/pdf/2005.03237.pdf

dfloat logAverage(dfloat a1, dfloat a2, dfloat loga1, dfloat loga2){

  dfloat da = a2-a1;
  dfloat aavg = .5*(a2+a1);
  dfloat f = da/aavg;
  dfloat v = f*f;

  if(fabs(f)<1e-4)
    return aavg*(1 + v*(-.2-v*(.0512 - v*0.026038857142857)));
  
  //  return -da/(log(a1)-log(a2));
  //  return -da/(log(a1/a2));
  return -da/(loga1-loga2); 
  
}


void entropy2conserved(dfloat gamma, dfloat *evar, dfloat *cvar){

  dfloat v1 = evar[0], v2 = evar[1], v3 = evar[2], v4 = evar[3];

  dfloat s  = p_gamma-v1 + 0.5*guardDivide(v2*v2+v3*v3,v4);  // OK TW HACK 2

  dfloat re = guardPow( (p_gamma-1.)/guardPow(-v4, p_gamma), 1./(p_gamma-1))*exp(-s/(p_gamma-1));  

  re = fabs(re); // TW: THIS IS SUPER CRITICAL
  
  dfloat r  = -re*v4; // OK
  dfloat ru =  re*v2; // OK
  dfloat rv =  re*v3; // OK
  dfloat E  =  re*(1- 0.5*guardDivide(v2*v2+v3*v3,v4));  
  
  cvar[0] = r;
  cvar[1] = ru;
  cvar[2] = rv;
  cvar[3] = E;
  
#if PRINT_WARNINGS==1
  if(re<=0 || v4>=0 || r<=0){
    printf("got small T=(re)=%g (evar=%g,%g,%g,%g) (s=%g) (qc=%g,%g,%g,%g)\n",
		      re, v1, v2, v3, v4, s, r, ru, rv, E);
  }
#endif

}


void conserved2entropy(dfloat gamma, dfloat *cvar, dfloat *evar){

  dfloat r = cvar[0], ru = cvar[1], rv = cvar[2], E = cvar[3];

  // TW HACK
  r = guardPositive(r);
  
#if PRINT_WARNINGS==1
  if(fabs(r)<1e-10) printf("r = %g, ru = %g, rv = %g, E = %g\n", r, ru, rv, E);
#endif
  
  dfloat u = guardDivide(ru,r);
  dfloat v = guardDivide(rv,r);

  dfloat re = E - 0.5*r*(u*u+v*v); // OK

  re = guardPositive(re);
  
  dfloat p = (p_gamma-1)*re;

  //  dfloat s = log(p/pow(r, gamma));  // OK
  dfloat s = guardLog( guardDivide(guardPositive(p), guardPow(r, p_gamma)) );
  //  dfloat s = log(p) - p_gamma*log(r);
		     
#if PRINT_WARNINGS==1
  if(p<1e-10) printf("#08: finding p=%g\n", p);
  if(r<1e-10) printf("#09: finding r=%g\n", r);
  if(re<1e-10) printf("#10: finding re=%g\n", re);
#endif
#if 1
  dfloat v1 = guardDivide(re*(p_gamma+1.-s)-E,re);   
  dfloat v2 = guardDivide(ru,re); 
  dfloat v3 = guardDivide(rv,re); 
  dfloat v4 = -guardPositive(guardDivide(r,re));  
#else
  dfloat reInv = 1./re;
  dfloat v1 = (p_gamma+1.-s)-E*reInv;   
  dfloat v2 = ru*reInv; 
  dfloat v3 = rv*reInv; 
  dfloat v4 = -r*reInv; 
#endif

  //  v4 = -fabs(v4); // TW
  // TW HACK
  v4 = -guardPositive(v4);
  
  evar[0] = v1;
  evar[1] = v2;
  evar[2] = v3;
  evar[3] = v4;
}

void conserved2primitive(dfloat gamma, dfloat *cvar, dfloat *pvar){

  dfloat r = cvar[0];

  dfloat rInv = guardDivide(1.,r);
  dfloat u = cvar[1]*rInv;
  dfloat v = cvar[2]*rInv;
  dfloat p = (p_gamma-1)*(cvar[3]-0.5*rInv*(cvar[1]*cvar[1]+cvar[2]*cvar[2]));

  pvar[0] = r;
  pvar[1] = u;
  pvar[2] = v;
  pvar[3] = p;
  
}

void primitive2conserved(dfloat gamma, dfloat *pvar, dfloat *cvar){
  cvar[0] = pvar[0]; // dens
  cvar[1] = pvar[0]*pvar[1]; // x-mom
  cvar[2] = pvar[0]*pvar[2]; // yx-mom
  cvar[3] = pvar[3]/(p_gamma-1.) + 0.5*pvar[0]*(pvar[1]*pvar[1]+pvar[2]*pvar[2]);
  
}

void entropy2primitive(dfloat gamma, dfloat *evar, dfloat *pvar){
  dfloat cvar[p_Nfields];
  entropy2conserved(gamma, evar, cvar);
  conserved2primitive(gamma, cvar, pvar);
}

void entropyStableTwoPointFlux(dfloat gamma, dfloat *evar1, dfloat *evar2,  dfloat *pvar1, dfloat *pvar2, dfloat *fx, dfloat *fy){

  // note the extra stack
  dfloat r1 = pvar1[0], u1 = pvar1[1], v1 = pvar1[2], p1 = pvar1[3], logr1 = pvar1[4], logbeta1 = pvar1[5], beta1 = pvar1[6];
  dfloat r2 = pvar2[0], u2 = pvar2[1], v2 = pvar2[2], p2 = pvar2[3], logr2 = pvar2[4], logbeta2 = pvar2[5], beta2 = pvar2[6];

  //  dfloat beta1 = 0.5*r1/p1; // OK
  //  dfloat beta2 = 0.5*r2/p2; // OK

  dfloat   rA = 0.5*(r1+r2); // OK
  dfloat   uA = 0.5*(u1+u2); // OK
  dfloat   vA = 0.5*(v1+v2); // OK

  dfloat  vnA = u1*u2+v1*v2;
  
  //  dfloat betaA    = 0.5*(beta1+beta2); // OK
  
  dfloat rAlog    =  logAverage(r1,    r2,    logr1,    logr2);    // OK
  dfloat betaAlog =  logAverage(beta1, beta2, logbeta1, logbeta2); // OK

  dfloat pA = guardDivide(rA,(beta1+beta2));  // OK

  dfloat EA = 0.5*guardDivide(rAlog,betaAlog*(p_gamma-1)) + 0.5*rAlog*vnA; // OK
  
  dfloat fx1 = rAlog*uA; // r*u OK
  dfloat fy1 = rAlog*vA; // r*v OK
  
  dfloat fx2 = rAlog*uA*uA + pA; //  r*u*u+p OK 
  dfloat fy2 = rAlog*vA*uA;      //  r*v*u   OK 
  
  dfloat fx3 = rAlog*uA*vA;      // r*u*v OK
  dfloat fy3 = rAlog*vA*vA + pA; // r*v*v+p OK

  // (logavg(r1,r2)/(2*(gamma-1)*logavg(0.5*r1/p1,0.5*r2/p2)  + pA + 0.5*r*(u*u+v*v)))
  dfloat fx4 = (EA + pA)*uA;   // OK
  dfloat fy4 = (EA + pA)*vA;   // OK

#if 0
  //
  fx1 += -lambda*(evar1[0]-evar2[0]);
  fx2 += -lambda*(evar1[1]-evar2[1]);
  fx3 += -lambda*(evar1[2]-evar2[2]);
  fx4 += -lambda*(evar1[3]-evar2[3]);

  fy1 += -lambda*(evar1[0]-evar2[0]);
  fy2 += -lambda*(evar1[1]-evar2[1]);
  fy3 += -lambda*(evar1[2]-evar2[2]);
  fy4 += -lambda*(evar1[3]-evar2[3]);
#endif
  
  fx[0] = fx1; fx[1] = fx2; fx[2] = fx3; fx[3] = fx4;
  fy[0] = fy1; fy[1] = fy2; fy[2] = fy3; fy[3] = fy4;
}

void averageTwoPointFlux(dfloat gamma, dfloat *pvar1, dfloat *pvar2, dfloat *fx, dfloat *fy){

  // note the extra stack
  dfloat r1 = pvar1[0], u1 = pvar1[1], v1 = pvar1[2], p1 = pvar1[3];
  dfloat r2 = pvar2[0], u2 = pvar2[1], v2 = pvar2[2], p2 = pvar2[3];

  dfloat E1 = p1/(p_gamma-1.) + 0.5*r1*(u1*u1+v1*v1);
  dfloat E2 = p2/(p_gamma-1.) + 0.5*r2*(u2*u2+v2*v2);
  
  fx[0] = 0.5*(r1*u1+r2*u2);
  fy[0] = 0.5*(r1*v1+r2*v2);

  fx[1] = 0.5*(r1*u1*u1+p1 +r2*u2*u2+p2);
  fy[1] = 0.5*(r1*v1*u1    +r2*v2*u2);

  fx[2] = 0.5*(r1*u1*v1    + r2*u2*v2);
  fy[2] = 0.5*(r1*v1*v1+p1 + r2*v2*v2+p2);

  fx[3] = 0.5*(u1*(E1+p1)  + u2*(E2+p2));
  fy[3] = 0.5*(v1*(E1+p1)  + v2*(E2+p2));
  
}

#if 1
void onePointFlux(dfloat gamma, dfloat *pvar1, dfloat *fx, dfloat *fy){

  dfloat r1 = pvar1[0], u1 = pvar1[1], v1 = pvar1[2], p1 = pvar1[3];

  dfloat E1 = p1/(p_gamma-1.) + 0.5*r1*(u1*u1+v1*v1);
  
  fx[0] = r1*u1;
  fy[0] = r1*v1;

  fx[1] = r1*u1*u1+p1;
  fy[1] = r1*v1*u1;

  fx[2] = r1*u1*v1;
  fy[2] = r1*v1*v1+p1;

  fx[3] = u1*(E1+p1);
  fy[3] = v1*(E1+p1);
  
}
#else
void onePointFlux(dfloat gamma, dfloat *cvar, dfloat *fx, dfloat *fy){
  
  dfloat r = cvar[0], ru = cvar[1], rv = cvar[2], E = cvar[3];
  
  dfloat rinv = 1./r;
  dfloat u = rinv*ru;
  dfloat v = rinv*rv;
  dfloat p = (p_gamma-1.)*(E - 0.5*r*(u*u+v*v));
  
  fx[0] = ru;
  fy[0] = rv;
  
  fx[1] = u*ru+p;
  fy[1] = v*ru;
  
  fx[2] = u*rv;
  fy[2] = v*rv+p;
  
  fx[3] = u*(E+p);
  fy[3] = v*(E+p);
}
#endif


#define ZERO 0.0
#define ONE  1.0
#define TWO  2.0
#define HALF 0.5

// https://home.cscamm.umd.edu/people/faculty/tadmor/pubs/files/Madrane_Tadmor_HYP2008.pdf
dfloat QH(dfloat a, dfloat deltaS){

  if(fabs(a)>0.5*deltaS) return fabs(a);

  return (a*a/deltaS) + 0.25*deltaS;
}

void esdgEntropyStableMatrixPenaltyTri2D(dfloat gamma, dfloat nx, dfloat ny, dfloat *cvar1, dfloat *cvar2, dfloat *penalty){

  dfloat var1[p_Nfields], var2[p_Nfields];
  
  // matrix dissipation flux
  conserved2entropy(p_gamma, cvar1, var1);
  conserved2entropy(p_gamma, cvar2, var2);
  
  const dfloat dv1 = var2[0]-var1[0];
  const dfloat dv2 = var2[1]-var1[1];
  const dfloat dv3 = var2[2]-var1[2];
  const dfloat dv4 = var2[3]-var1[3];	

  conserved2primitive(p_gamma, cvar1, var1);
  conserved2primitive(p_gamma, cvar2, var2);

  const dfloat r1 = fabs(var1[0]), r2 = fabs(var2[0]);
  const dfloat u1 = var1[1],       u2 = var2[1];
  const dfloat v1 = var1[2],       v2 = var2[2];
  const dfloat p1 = fabs(var1[3]), p2 = fabs(var2[3]);

  const dfloat beta1 = fabs(0.5*r1/p1);
  const dfloat beta2 = fabs(0.5*r2/p2);
  
  const dfloat u2avg = 0.5*(u1*u1+u2*u2);
  const dfloat v2avg = 0.5*(v1*v1+v2*v2);
  const dfloat uavg  = 0.5*(u1+u2);
  const dfloat vavg  = 0.5*(v1+v2);
  const dfloat betalog = fabs(logAverage(beta1,beta2,log(beta1),log(beta2)));
  const dfloat pavg    = 0.5*(r1+r2)/(TWO*0.5*(beta1+beta2));
  const dfloat rholog  = fabs(logAverage(r1,r2,log(r1),log(r2)));
  
  const dfloat ubar = TWO*(uavg*uavg+vavg*vavg)-(u2avg+v2avg);
  const dfloat un   = (uavg*nx+vavg*ny);
  const dfloat h    = p_gamma/(TWO*betalog*(p_gamma-ONE)) + HALF*ubar;
  const dfloat a    = sqrt(p_gamma*pavg/rholog);
  
  const dfloat R11 = ONE;
  const dfloat R12 = ONE;
  const dfloat R13 = ZERO;
  const dfloat R14 = ONE;	
  
  const dfloat R21 = uavg-a*nx;
  const dfloat R22 = uavg;
  const dfloat R23 = ny;
  const dfloat R24 = uavg+a*nx;	
  
  const dfloat R31 = vavg-a*ny;
  const dfloat R32 = vavg;
  const dfloat R33 = -nx;
  const dfloat R34 = vavg+a*ny;
  
  const dfloat R41 = h-a*un;
  const dfloat R42 = HALF*ubar;
  const dfloat R43 = uavg*ny - vavg*nx;
  const dfloat R44 = h+a*un;

#if 1
  const dfloat D11 = fabs(un - a)*rholog/(TWO*p_gamma);
  const dfloat D22 = fabs(un)*rholog*(p_gamma-ONE)/p_gamma;
  const dfloat D33 = fabs(un)*pavg;
  const dfloat D44 = fabs(un + a)*rholog/(TWO*p_gamma);
#else
  // Harten correction to control nearly stationnary waves
  dfloat delt = .25;
  dfloat deltaS = delt*(fabs(uavg)+fabs(vavg)+a);

  const dfloat D11 = QH(un - a, deltaS)*rholog/(TWO*p_gamma);
  const dfloat D22 = QH(un,     deltaS)*rholog*(p_gamma-ONE)/p_gamma;
  const dfloat D33 = QH(un,     deltaS)*pavg;
  const dfloat D44 = QH(un + a, deltaS)*rholog/(TWO*p_gamma);
#endif
  // R*D*R' * [[v]]
  const dfloat f1 = D11*(R11*dv1 + R21*dv2 + R31*dv3 + R41*dv4);
  const dfloat f2 = D22*(R12*dv1 + R22*dv2 + R32*dv3 + R42*dv4);
  const dfloat f3 = D33*(R13*dv1 + R23*dv2 + R33*dv3 + R43*dv4);
  const dfloat f4 = D44*(R14*dv1 + R24*dv2 + R34*dv3 + R44*dv4);	
  
  penalty[0] = R11*f1 + R12*f2 + R13*f3 + R14*f4;
  penalty[1] = R21*f1 + R22*f2 + R23*f3 + R24*f4;
  penalty[2] = R31*f1 + R32*f2 + R33*f3 + R34*f4;
  penalty[3] = R41*f1 + R42*f2 + R43*f3 + R44*f4;

#if 0 
  penalty[0] /= PSCALE;
  penalty[1] /= PSCALE;
  penalty[2] /= PSCALE;
  penalty[3] /= PSCALE;
#endif
}

void esdgEntropyStableJacobianMatrixPenaltyTri2D(dfloat gamma,
						  dfloat nx, dfloat ny, 
						  dfloat *pvar1,
						  dfloat *pvar2,
						  dfloat *evar1,
						  dfloat *evar2,
						  dfloat *penalty){

  // https://liu.diva-portal.org/smash/get/diva2:1315778/FULLTEXT01
  
  dfloat var1[p_Nfields], var2[p_Nfields];

  dfloat r1 = fabs(pvar1[0]), u1 = pvar1[1], v1 = pvar1[2], p1 = fabs(pvar1[3]);
  dfloat r2 = fabs(pvar2[0]), u2 = pvar2[1], v2 = pvar2[2], p2 = fabs(pvar2[3]);

  dfloat dv1 = evar2[0] - evar1[0];
  dfloat dv2 = evar2[1] - evar1[1];
  dfloat dv3 = evar2[2] - evar1[2];
  dfloat dv4 = evar2[3] - evar1[3];
  
  dfloat rLOG = logAverage(r1, r2, log(r1), log(r2));
  dfloat uA = 0.5*(u1+u2);
  dfloat vA = 0.5*(v1+v2);

  dfloat u2A = 0.5*(u1*u1+u2*u2);
  dfloat v2A = 0.5*(v1*v1+v2*v2);
  
  dfloat pA = 0.5*(p1+p2);

  dfloat beta1 = r1/(2.*p1);
  dfloat beta2 = r2/(2.*p2);
  dfloat betaLOG = logAverage(beta1, beta2, log(beta1), log(beta2));
  
  dfloat pLOG = rLOG/(2.*betaLOG);
  dfloat u2BAR = 2*(uA*uA+vA*vA)-u2A-v2A;
  
  dfloat EBAR = pLOG/(p_gamma-1) + 0.5*rLOG*u2BAR;
  
  dfloat a1 = sqrt(abs(p_gamma*p1/r1));
  dfloat a2 = sqrt(abs(p_gamma*p2/r2));

  dfloat un1 = nx*u1+ny*v1;
  dfloat un2 = nx*u2+ny*v2;
  dfloat unA = 0.5*(un1+un2);
#if 0	 
  
  dfloat lambda1 = abs(un1) + a1;
  dfloat lambda2 = abs(un2) + a2;
  
  dfloat lambdaMAX = (lambda1>lambda2) ? lambda1:lambda2;
#else
  dfloat rA =  0.5*(r1+r2);
  const dfloat aA = sqrt(fabs(p_gamma*pA/rA));
  dfloat lambdaMAX = fabs(unA) + aA;
#endif

  dfloat Hhat44 = ((pLOG*pLOG)/(p_gamma-1) + EBAR*EBAR)/rLOG + pA*(uA*uA+vA*vA);
  
  penalty[0] =    rLOG*dv1 +         rLOG*uA*dv2 +         rLOG*vA*dv3 +         EBAR*dv4;
  penalty[1] = rLOG*uA*dv1 + (rLOG*uA*uA+pA)*dv2 + (rLOG*uA*vA   )*dv3 + uA*(EBAR+pA)*dv4;
  penalty[2] = rLOG*vA*dv1 + (rLOG*vA*uA   )*dv2 + (rLOG*vA*vA+pA)*dv3 + vA*(EBAR+pA)*dv4;
  penalty[3] =    EBAR*dv1 + uA*(EBAR+pA   )*dv2 +    vA*(EBAR+pA)*dv3 +       Hhat44*dv4;

  penalty[0] *= lambdaMAX;
  penalty[1] *= lambdaMAX;
  penalty[2] *= lambdaMAX;
  penalty[3] *= lambdaMAX;
}
