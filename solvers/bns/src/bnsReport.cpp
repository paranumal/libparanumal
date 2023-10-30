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


#define TGV 0


#include "bns.hpp"

void bns_t::Report(dfloat time, int tstep){

  static int frame=0;

  //compute q.M*q
  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
  deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
  mesh.MassMatrixApply(o_q, o_Mq);

  dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));
  o_Mq.free();

  if(mesh.rank==0)
    printf("%5.2f (%d), %5.2f (time, timestep, norm)\n", time, tstep, norm2);

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

    //compute vorticity
    deviceMemory<dfloat> o_Vort = platform.reserve<dfloat>(mesh.dim*mesh.Nelements*mesh.Np);
    vorticityKernel(mesh.Nelements, mesh.o_vgeo, mesh.o_D, o_q, c, o_Vort);

    memory<dfloat> Vort(mesh.dim*mesh.Nelements*mesh.Np);

    // copy data back to host
    o_q.copyTo(q);
    o_Vort.copyTo(Vort);

    // output field files
    std::string name;
    settings.getSetting("OUTPUT FILE NAME", name);
    char fname[BUFSIZ];
    sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

    PlotFields(q, Vort, std::string(fname));
  }

  o_q.copyTo(q); ComputeForces(time); 


  // #endif
}




void bns_t::ComputeForces(const dfloat T){

  // Hard coded weights i.e. sum(MM1D,1), MM1D = inv(V1)'*V1
  
  memory<dfloat> W; W.malloc(mesh.Nfp, 0);

  memory<dfloat> Mf; Mf.malloc(mesh.Nfp*mesh.Nfp,0); 

  if(mesh.dim==2){
    const int Np = (mesh.N+1); 
    memory<dfloat> r, V;
    mesh.Nodes1D(mesh.N, r); //Gauss-Legendre-Lobatto nodes
    mesh.Vandermonde1D(mesh.N, r, V); 
    mesh.MassMatrix1D(Np, V, Mf);   
  }else if(mesh.dim==3){
    const int Np = (mesh.N+1)*(mesh.N+2)/2; 
    memory<dfloat> r, s, V;
    mesh.NodesTri2D(mesh.N, r, s); //Gauss-Legendre-Lobatto nodes
    mesh.VandermondeTri2D(mesh.N, r, s, V); 
    mesh.MassMatrixTri2D(Np, V, Mf);   
  }

  for(int i=0; i< mesh.Nfp; i++){
    dfloat sum = 0; 
      for(int j=0; j< mesh.Nfp; j++){
        sum += Mf[i + j*mesh.Nfp]; 
      }
      W[i] = sum; 
      // printf(" %.8f \n", W[i]);
  }
  //
  dfloat FVx = 0.0, FVy = 0.0, FVz = 0.0; 
  dfloat FPx = 0.0, FPy = 0.0, FPz = 0.0; 
   for(int e=0;e<mesh.Nelements;++e){
  
    for(int f=0;f<mesh.Nfaces;++f){
  
        int bc = mesh.EToB[e*mesh.Nfaces+f];
  
        if(bc==1){ // HarCoded here i.e. wall
  
          for(int n=0;n<mesh.Nfp; n++){
            dfloat nx=0.0, ny=0.0, nz=0.0, sJ = 0.0; 
            const dlong sid = e*mesh.Nfaces*mesh.Nsgeo + f*mesh.Nsgeo;
            nx = mesh.sgeo[sid+mesh.NXID];
            ny = mesh.sgeo[sid+mesh.NYID];

            nz = mesh.dim==3 ? mesh.sgeo[sid+mesh.NZID]: 0.0;

            sJ = mesh.sgeo[sid+mesh.SJID];

            const dlong vid  = e*mesh.Nfp*mesh.Nfaces + f*mesh.Nfp + n;
            const dlong idM  = mesh.vmapM[vid];
            // load traces
            const int vidM = idM%mesh.Np;

            const dlong qidM = e*mesh.Np*Nfields + vidM;

            // Read trace values
            dfloat q1  = 0.0, q2  = 0.0, q3  = 0.0, q4  = 0.0, q5  = 0.0; 
            dfloat q6  = 0.0, q7  = 0.0, q8  = 0.0, q9  = 0.0, q10 = 0.0; 

            // Read trace values
            q1  = q[qidM + 0*mesh.Np];
            q2  = q[qidM + 1*mesh.Np];
            q3  = q[qidM + 2*mesh.Np];
            q4  = q[qidM + 3*mesh.Np];
            q5  = q[qidM + 4*mesh.Np];
            q6  = q[qidM + 5*mesh.Np];

            if(mesh.dim==3){
              q7  = q[qidM + 6*mesh.Np];
              q8  = q[qidM + 7*mesh.Np];
              q9  = q[qidM + 8*mesh.Np];
              q10 = q[qidM + 9*mesh.Np];
           }


           dfloat P, s11, s12, s13, s22, s23, s33; 

           if(mesh.dim==2){
            P  =   RT*q1; 
            s11 = -RT*(sqrt(2) * q5 - q2*q2/ q1);
            s22 = -RT*(sqrt(2) * q6 - q3*q3/ q1);
            s12 = -RT*(          q4 - q2*q3/ q1);
            
            FVx += -W[n]*sJ*(s11*nx + s12*ny);   
            FVy += -W[n]*sJ*(s12*nx + s22*ny);   
            FPx +=  W[n]*sJ*(P*nx);   
            FPy +=  W[n]*sJ*(P*ny);  

           }else if(mesh.dim==3){

            P  =   RT*q1; 

            s11 = -RT*(sqrt(2) * q5 - q2*q2/ q1);
            s22 = -RT*(sqrt(2) * q6 - q3*q3/ q1);
            s33 = -RT*(sqrt(2) * q7 - q4*q4/ q1);
            
            s12 = -RT*(          q8 - q2*q3/ q1);
            s13 = -RT*(          q9 - q2*q4/ q1);
            s23 = -RT*(          q10- q3*q4/ q1);

            FVx += -W[n]*sJ*(s11*nx + s12*ny + s13*nz);   
            FVy += -W[n]*sJ*(s12*nx + s22*ny + s23*nz);   
            FVz += -W[n]*sJ*(s13*nx + s23*ny + s33*nz);   
            FPx +=  W[n]*sJ*(P*nx);   
            FPy +=  W[n]*sJ*(P*ny);  
            FPz +=  W[n]*sJ*(P*nz);  
           }           
          }
        }
      }
}

  dfloat gFVx = FVx; 
  dfloat gFVy = FVy; 
  dfloat gFPx = FPx; 
  dfloat gFPy = FPy; 

   comm.Allreduce(gFVx, Comm::Sum);
   comm.Allreduce(gFVy, Comm::Sum);
   comm.Allreduce(gFPx, Comm::Sum);
   comm.Allreduce(gFPy, Comm::Sum);

   dfloat gFVz = 0.0, gFPz = 0.0;  
if(mesh.dim==3){
    gFVz = FVz;  
    gFPz = FPz; 
    comm.Allreduce(gFVz, Comm::Sum);
    comm.Allreduce(gFPz, Comm::Sum);
}



  if(mesh.rank==0){
    char fname[BUFSIZ];
    sprintf(fname, "BNSForceData_N%d.dat", mesh.N);
    FILE *fp; fp = fopen(fname, "a");  
    if(mesh.dim==2)
      fprintf(fp, "%.4e %.8e %.8e %.8e %.8e\n", T, gFVx,  gFVy,  gFPx,  gFPy); 
    
    if(mesh.dim==3)
      fprintf(fp, "%.4e %.8e %.8e %.8e %.8e %.8e %.8e\n", T, gFVx,  gFVy, gFVz, gFPx, gFPy, gFPz); 

    fclose(fp);
  }



#if TGV==1

  // Compute Kinetic Energy

  memory<dfloat> U; U.malloc(mesh.Nelements*mesh.Np*Nfields, 0.0); 

  for(int e=0; e< mesh.Nelements; e++){
    for(int n=0; n<mesh.Np; n++){
      const dlong id = e*mesh.Np*Nfields + n; 
      
      dfloat q1 = 0.0, q2= 0.0, q3= 0.0, q4= 0.0; 

      q1 = q[id + 0*mesh.Np]; 
      q2 = q[id + 1*mesh.Np]; 
      q3 = q[id + 2*mesh.Np];

      q4 = (mesh.dim==3) ? q[id + 3*mesh.Np]: 0.0;

      const dfloat u   = q2/q1; 
      const dfloat v   = q3/q1; 
      const dfloat w   = q4/q1; 
      //
      U[e*mesh.Np*Nfields + n + 0*mesh.Np] =  q2/q1; 
      U[e*mesh.Np*Nfields + n + 1*mesh.Np] =  q3/q1; 
      U[e*mesh.Np*Nfields + n + 2*mesh.Np] = (mesh.dim==3) ? q4/q1 : 0.0;
    }
  }


  deviceMemory<dfloat> o_U; o_U = platform.malloc<dfloat>(U);  
//compute q.M*q
  mesh.MassMatrixApply(o_U, o_Mq);

  // dfloat iV = (mesh.dim==3) ? 1.0/pow(2.0*M_PI, 3): 1.0/pow(2.0*M_PI, 2);
  dfloat iV = 1.0;

  dlong Nentries = mesh.Nelements*mesh.Np*Nfields;

  dfloat ke  = 0.5*iV*platform.linAlg().innerProd(Nentries, o_U, o_Mq, mesh.comm);

  printf("ke = %.8e \n", ke);

  if(mesh.rank==0){
    char fname[BUFSIZ];
    sprintf(fname, "TGV_Energy_N%d.dat", mesh.N);
    FILE *fp; fp = fopen(fname, "a");  
    
    fprintf(fp, "%.4e %.8e \n", T, ke); 

    fclose(fp);
  }


#endif







  }






