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

#include "esdg.hpp"

void esdg_t::Report(dfloat time, int tstep){

  static int frame=0;

  //compute vorticity
  
  //compute q.M*q
  MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_q, o_Mq);

  // check o_LIFTT => o_LIFT substitution
  // check o_Dmatrices => o_D
  // TW fix this later ?
  if(1) 
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D,
		    o_q, o_Vort);
  else
    dgVorticityKernel(mesh.Nelements, mesh.o_vmapM, mesh.o_vmapP, mesh.o_vgeo, mesh.o_sgeo, mesh.o_D, mesh.o_LIFT, 
		    o_q, o_Vort);

  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  dfloat norm2 = 0;

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
    o_q.copyTo(q);
    o_Vort.copyTo(Vort);
    
    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_E%06d_N%02d_%04d_%04d.vtu", name.c_str(), mesh.Nelements, mesh.N, mesh.rank, frame++);
    
    PlotFields(q, Vort, fname);
  }
  
  
  if(1){
    
    // copy data back to host
    o_q.copyTo(q);
    // compute mach number
    dfloat minMach = 1e9, maxMach = 0;
    dfloat minDens = 1e9, maxDens = 0;
    dfloat minPres = 1e9, maxPres = 0;

    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
	dlong base = e*mesh.Np*Nfields+n;
	
	dfloat r  = q[base+0*mesh.Np];
	dfloat ru = q[base+1*mesh.Np];
	dfloat rv = q[base+2*mesh.Np];
	dfloat rw = 0, EN = 0;
	if(mesh.dim==2){
	  rw = 0;
	  EN = q[base+3*mesh.Np];
	}else{
	  rw = q[base+3*mesh.Np];
	  EN = q[base+4*mesh.Np];
	}

	dfloat p = (gamma-1)*(EN-0.5*(ru*ru+rv*rv+rw*rw)/r);
	dfloat a = sqrt(gamma*p/r);
	dfloat magU = sqrt(ru*ru+rv*rv+rw*rw)/r;
	dfloat mach = magU/a;
	
	minMach = mymin(minMach, mach);
	maxMach = mymax(maxMach, mach);
	minDens = mymin(minDens, r);
	maxDens = mymax(maxDens, r);
	minPres = mymin(minPres, p);
	maxPres = mymax(maxPres, p);

      }
    }

    printf("MACH NUMBER in RANGE [%7.5E,%7.5E]\n", minMach, maxMach);
    printf("DENSITY     in RANGE [%7.5E,%7.5E]\n", minDens, maxDens);
    printf("PRESSURE    in RANGE [%7.5E,%7.5E]\n", minPres, maxPres);

    int orderPlus = mesh.cubN;

    printf("%d %d % 7.5E % 7.5E % 7.5E % 7.5E % 7.5E % 7.5E % 7.5E",
	   mesh.N, orderPlus, time,  minMach, maxMach, minDens, maxDens, minPres, maxPres);

    if(cubature)
      printf("%%%% degree, cubatureN, filter recon degree, time, minMach, maxMach, minDens, maxDens, minPres, maxPres \n");
    else
      printf("%%%% degree, degree_SN, flux degree, time, minMach, maxMach, minDens, maxDens, minPres, maxPres \n");
    
    fflush(stdout);
  }

  if(1){

    // TW - check contents of o_cubInterp (was o_cubInterpT)
    errorKernel(mesh.Nelements, gamma, time,
		mesh.o_vgeo, o_cx, o_cy, o_cz,
		mesh.o_cubx, mesh.o_cuby, mesh.o_cubz, mesh.o_cubInterp,
		o_cubw, o_q, o_linfError, o_l2Error);
    
    o_linfError.copyTo(linfError);
    o_l2Error.copyTo(l2Error);

    memory<dfloat> l2ErrorTotal, linfErrorTotal;
    l2ErrorTotal.malloc(Nfields);
    linfErrorTotal.malloc(Nfields);

    for(dlong e=0;e<mesh.Nelements;++e){
      for(int f=0;f<Nfields;++f){
	l2ErrorTotal[f]   += l2Error[e*Nfields+f];
	linfErrorTotal[f]  = mymax(linfErrorTotal[f], linfError[e*Nfields+f]);
      }
    }

    // need MPI reduce here

    int orderPlus =  mesh.cubN;
    
    printf("%d %d %d %7.5E ", mesh.Nelements, mesh.N, orderPlus, time);

    // print out LINF errors
    for(int f=0;f<Nfields;++f){
      printf(" % 7.5E ", linfErrorTotal[f]);
    }

    // print out L2 errors
    for(int f=0;f<Nfields;++f){
      printf(" % 7.5E ", sqrt(l2ErrorTotal[f]));
    }

    if(cubature)
      printf("%%%% Nelements, degree, cubatureN, flux degree, time, MAX ERRORSx4, L2 ERRORSx4   \n");
    else
      printf("%%%% Nelements, degree, degree_SN, flux degree, time, MAX ERRORSx4, L2 ERRORSx4   \n");
      
    fflush(stdout);
  }

  if (settings.compareSetting("CHECKPOINT SAVE","TRUE")) {
    o_q.copyTo(q);
    saveCheckpoint(q, time);
  }
}
