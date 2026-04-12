/*

The MIT License (MIT)

Copyright (c) 2017-2026 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "cns.hpp"
void cns_t::setupProbe(){

	reportProbes = settings.compareSetting("REPORT PROBES", "TRUE") ? 1:0;
	if (!reportProbes) return;

// Read probe data
	readProbe(); 

// Locate elements and find local coordinates
	switch (mesh.elementType) {
	case Mesh::TRIANGLES:
		locateProbesTri2D();
		break;
	case Mesh::QUADRILATERALS:
		locateProbesQuad2D();
		break;
	case Mesh::TETRAHEDRA:
		locateProbesTet3D();
		break;
	case Mesh::HEXAHEDRA:
		locateProbesHex3D();
		break;
	}

	int NprobeLocated = 0;
	// Find the number of unlocated probes
	mesh.comm.Allreduce(NprobeLocal, NprobeLocated, MPI_SUM);
	mesh.comm.Allreduce(probeB,MPI_SUM);

	if(NprobeLocated!=NprobeGlobal){
		if(mesh.rank==0){
			printf("Initially number of probes located %d over %d total probes\n", NprobeLocated, NprobeGlobal);
			printf("Locating missing probes--------- :");
		}
		memory<dfloat> bestDG(NprobeGlobal); 
		mesh.comm.Allreduce(bestD,bestDG, MPI_MIN);
		for(int p=0; p<NprobeGlobal; p++){
			if(probeB[p]==0 && abs(bestD[p] - bestDG[p])<1e-12){
				NprobeLocal++; 
				for(int d=0; d<mesh.dim;d++){
				  probeR[NprobeLocal*mesh.dim+d] = bestR[p*mesh.dim+d]; 
				}
				probeE[NprobeLocal] = bestE[p]; 
				probeIDl[NprobeLocal] = probeIDa[p]; 
			}
		}

		// Find the number of unlocated probes
	mesh.comm.Allreduce(NprobeLocal, NprobeLocated, MPI_SUM);

	if(mesh.rank==0){
			printf("done\n");
		}
	}


  if(mesh.rank==0){
		printf("Number of probes located %d over %d total probes\n", NprobeLocated, NprobeGlobal);
	}

	probeR.realloc(NprobeLocal*mesh.dim); 
	probeE.realloc(NprobeLocal); 
	probeIDl.realloc(NprobeLocal); 

	// for(int p=0; p<NprobeLocal;p++){
	// 		printf("rank=%d %d %.4f %.4f\n", mesh.rank, p, probeR[p*mesh.dim+0],probeR[p*mesh.dim+1]);
	// }

// Build Interpolation Matrix  
	switch (mesh.elementType) {
	case Mesh::TRIANGLES:
		interpolateProbesTri2D();
		break;
	case Mesh::QUADRILATERALS:
		interpolateProbesQuad2D();
		break;
	case Mesh::TETRAHEDRA:
		interpolateProbesTet3D();
		break;
	case Mesh::HEXAHEDRA:
		interpolateProbesHex3D();
		break;
	}

// MPI setup
	if(mesh.rank == 0){
		probeRecvCount.malloc(mesh.size);
		probeRecvOffset.malloc(mesh.size);
	}

// Collect local number of probes in every processor
	memory<int> NprobeL(1); NprobeL[0] = NprobeLocal;
	mesh.comm.Gather(NprobeL, probeRecvCount, 0, 1);
//
	if(mesh.rank == 0){
		probeRecvOffset[0] = 0;
		for(int i=0;i<mesh.size;++i){
			if(i>0) 
				probeRecvOffset[i] = probeRecvOffset[i-1]+probeRecvCount[i-1];
		}
	}


}







void cns_t::readProbe(){
	std::string probeInputFile; 
	settings.getSetting("PROBE INPUT FILE",      probeInputFile);
	std::ifstream fp(probeInputFile);
	if (!fp.is_open()) {
		printf("Could not open probe file: %s\n", probeInputFile.c_str());
		LIBP_ABORT("ReadProbeile: could not open probe file", 1);
	}

	std::string line;
	NprobeGlobal =0;
	while (std::getline(fp, line)) {
		if (line.empty() || line[0] == '#') continue;

		std::istringstream iss(line);
		dfloat x, y, z;
		if (mesh.dim==3) {
			if (iss >> x >> y >> z) ++NprobeGlobal;
		} else {
			if (iss >> x >> y) ++NprobeGlobal;
		}
	}
	fp.close();
	probeX.malloc(NprobeGlobal*mesh.dim);
	// All probe IDs
	probeIDa.malloc(NprobeGlobal);



	fp.open(probeInputFile);
	if (!fp.is_open()) {
		printf("could not reopen probe file: %s\n", probeInputFile.c_str());
		LIBP_ABORT("ReadProbeile: could not reopen probe file", 1);
	}

	dlong n = 0;
	while (std::getline(fp, line)) {
		if (line.empty() || line[0] == '#') 
			continue;

		std::istringstream iss(line);
		dfloat x, y, z; int id; 
		if (mesh.dim==3){ 
			if (!(iss >> id>>x >> y >> z)) continue;
		} 
		else { 
			if (!(iss >> id>> x >> y)) continue;
		}
		probeX[n*mesh.dim+0] = x;
		probeX[n*mesh.dim+1] = y;
		if(mesh.dim==3) { 
			probeX[n*mesh.dim+2] = z;
		}
		probeIDa[n]  = id;
		++n;
	}
	fp.close();

	// if(mesh.rank==0){
	// 	for(n=0; n<NprobeGlobal; n++){
	// // printf("%d, %.4f %.4f %.4f\n", probeID[n], probeX[n*mesh.dim + 0], probeX[n*mesh.dim + 1], probeX[n*mesh.dim + 2]);
	// 		printf("%d, %.4f %.4f\n", probeIDa[n], probeX[n*mesh.dim + 0], probeX[n*mesh.dim + 1]);
	// 	}
	// }

}



void cns_t::probeInterp(const memory<dfloat> u, memory<dfloat> Iu){
// interpolate
	for(int n=0;n<NprobeLocal;++n){
		dfloat un = 0;
		for(int m=0;m<mesh.Np;++m){
			un += probeI[n*mesh.Np+m]*u[m];
		}
		Iu[n] = un;
	}
}

void cns_t::interpolateProbesTri2D(){
	probeI.malloc(NprobeLocal*mesh.Np); 

	memory<dfloat> r(NprobeLocal); 
	memory<dfloat> s(NprobeLocal); 
	for(int p=0; p<NprobeLocal; p++){
		r[p] = probeR[p*mesh.dim +0];
		s[p] = probeR[p*mesh.dim +1];
	}
	mesh.InterpolationMatrixTri2D(mesh.N, mesh.r, mesh.s, r, s, probeI); 
}


void cns_t::interpolateProbesTet3D(){
	probeI.malloc(NprobeLocal*mesh.Np); 

	memory<dfloat> r(NprobeLocal); 
	memory<dfloat> s(NprobeLocal); 
	memory<dfloat> t(NprobeLocal); 

	for(int p=0; p<NprobeLocal; p++){
		r[p] = probeR[p*mesh.dim +0];
		s[p] = probeR[p*mesh.dim +1];
		t[p] = probeR[p*mesh.dim +2];
	}
	mesh.InterpolationMatrixTet3D(mesh.N, mesh.r, mesh.s,  mesh.s, r, s, t, probeI); 
}


void cns_t::interpolateProbesQuad2D(){

}


void cns_t::interpolateProbesHex3D(){
}



void cns_t::locateProbesTri2D(){
	probeR.malloc(NprobeGlobal*mesh.dim); 
	probeE.malloc(NprobeGlobal); 

	bestE.malloc(NprobeGlobal); 
	bestD.malloc(NprobeGlobal); 
	bestR.malloc(NprobeGlobal*mesh.dim); 

	probeIDl.malloc(NprobeGlobal);
	// Check that this probe is fixed
	probeB.calloc(NprobeGlobal); 

	memory<dfloat> A((mesh.dim+1)*mesh.Nverts ); 
	memory<dfloat> b((mesh.dim+1)*NprobeGlobal); 
	memory<dfloat> c((mesh.dim+1)*NprobeGlobal);


	NprobeLocal = 0; 
	dfloat tol = 1e-12; 

	for(int p=0; p<NprobeGlobal; p++){
		bestR[p*mesh.dim + 0] = 1e12; 
		bestR[p*mesh.dim + 1] = 1e12; 
		bestE[p] = 0; 

	}

	// fill up RHS i.e. Ac = b
	for(int p=0; p<NprobeGlobal; p++){
		b[p*(mesh.dim+1) + 0] = 1.0; 
		b[p*(mesh.dim+1) + 1] = probeX[p*mesh.dim+0]; 
		b[p*(mesh.dim+1) + 2] = probeX[p*mesh.dim+1]; 
	}

	for(dlong e=0; e<mesh.Nelements; e++){

		for (int v=0;v<mesh.Nverts;v++) {
			A[v*mesh.Nverts + 0] = 1.0;
			A[v*mesh.Nverts + 1] = mesh.EX[e*mesh.Nverts+v];
			A[v*mesh.Nverts + 2] = mesh.EY[e*mesh.Nverts+v];
		} 

		linAlg_t::matrixRightSolve(NprobeGlobal,(mesh.dim+1),  b, mesh.Nverts,(mesh.dim+1),  A,  c);

		for(int p=0; p<NprobeGlobal; p++){
			const int pid = probeIDa[p];

			if(probeB[p]==0){

				dfloat l1 = c[p*(mesh.dim+1) + 2]; 
				dfloat l2 = c[p*(mesh.dim+1) + 0]; 
				dfloat l3 = c[p*(mesh.dim+1) + 1]; 

				dfloat lmin = std::min(l1, std::min(l2,l3)); 

				const dfloat r = 2.0*l3-1.0; // r
				const dfloat s = 2.0*l1-1.0; // r

				if(lmin>tol){
					probeR[NprobeLocal*mesh.dim + 0] = r;  
			  	probeR[NprobeLocal*mesh.dim + 1] = s; // s
			  	probeE[NprobeLocal] = e;
			  	probeIDl[NprobeLocal] = pid; 
			   	probeB[p] = 1; // fix this prope
			   	NprobeLocal++; 
		   	break; 
		   }else{
					const dfloat d1 = (s < -1.0) ? (-1.0 - s) : 0.0;
  				const dfloat d2 = (r < -1.0) ? (-1.0 - r) : 0.0;
  				const dfloat d3 = (r + s > 0.0) ? (r + s) : 0.0;
  				const dfloat dist2 = 0.25*(d1*d1 + d2*d2 + d3*d3);

  				const dfloat rc =bestR[p*mesh.dim + 0]; 
  				const dfloat sc =bestR[p*mesh.dim + 1]; 

					const dfloat dc1 = (sc < -1.0) ? (-1.0 - sc) : 0.0;
  				const dfloat dc2 = (rc < -1.0) ? (-1.0 - rc) : 0.0;
  				const dfloat dc3 = (rc + sc > 0.0) ? (rc + sc) : 0.0;
  				const dfloat distc2 = 0.25*(dc1*dc1 + dc2*dc2 + dc3*dc3);

		   	  bestD[p] = dist2<distc2? dist2:distc2 ; // r
		   	  bestR[p*mesh.dim + 0] = dist2<distc2? r: rc ; // r
			  	bestR[p*mesh.dim + 1] = dist2<distc2? s: sc ; // r
			  	bestE[p]              = dist2<distc2? e: bestE[p];  
		   }
		 }

		}
	}

	

	// for(int p=0; p<NprobeLocal; p++){
	// 	printf("%d %d %d %.4e %.4e\n ",mesh.rank, probeE[p],  probeIDl[p], probeR[p*mesh.dim+0],probeR[p*mesh.dim+1]); 
	// }
}

void cns_t::locateProbesTet3D(){
	probeR.malloc(NprobeGlobal*mesh.dim); 
	probeE.malloc(NprobeGlobal); 

	bestE.malloc(NprobeGlobal); 
	bestD.malloc(NprobeGlobal); 
	bestR.malloc(NprobeGlobal*mesh.dim); 

	probeIDl.malloc(NprobeGlobal);
	// Check that this probe is fixed
	probeB.calloc(NprobeGlobal); 

	memory<dfloat> A((mesh.dim+1)*mesh.Nverts ); 
	memory<dfloat> b((mesh.dim+1)*NprobeGlobal); 
	memory<dfloat> c((mesh.dim+1)*NprobeGlobal);


	NprobeLocal = 0; 
	dfloat tol = 1e-12; 
	for(int p=0; p<NprobeGlobal; p++){
		bestR[p*mesh.dim + 0] = 1e12; 
		bestR[p*mesh.dim + 1] = 1e12; 
		bestR[p*mesh.dim + 2] = 1e12; 
		bestE[p] = 0; 

	}

	// fill up RHS i.e. Ac = b
	for(int p=0; p<NprobeGlobal; p++){
		b[p*(mesh.dim+1) + 0] = 1.0; 
		b[p*(mesh.dim+1) + 1] = probeX[p*mesh.dim+0]; 
		b[p*(mesh.dim+1) + 2] = probeX[p*mesh.dim+1]; 
		b[p*(mesh.dim+1) + 3] = probeX[p*mesh.dim+2]; 
	}

	for(dlong e=0; e<mesh.Nelements; e++){

		for (int v=0;v<mesh.Nverts;v++) {
			A[v*mesh.Nverts + 0] = 1.0;
			A[v*mesh.Nverts + 1] = mesh.EX[e*mesh.Nverts+v];
			A[v*mesh.Nverts + 2] = mesh.EY[e*mesh.Nverts+v];
			A[v*mesh.Nverts + 3] = mesh.EZ[e*mesh.Nverts+v];
		} 

		linAlg_t::matrixRightSolve(NprobeGlobal,(mesh.dim+1),  b, mesh.Nverts,(mesh.dim+1),  A,  c);

	
		for(int p=0; p<NprobeGlobal; p++){
			const int pid = probeIDa[p];
			if(probeB[p] ==0){
				dfloat l1 = c[p*(mesh.dim+1) + 3]; 
				dfloat l2 = c[p*(mesh.dim+1) + 2]; 
				dfloat l3 = c[p*(mesh.dim+1) + 0]; 
				dfloat l4 = c[p*(mesh.dim+1) + 1]; 

				dfloat lmin = std::min(l1,std::min(l2, std::min(l3,l4))); 

				const dfloat r = 2.0*l4-1.0;
				const dfloat s = 2.0*l2-1.0;
				const dfloat t = 2.0*l1-1.0;

				if(lmin>tol){
					probeR[NprobeLocal*mesh.dim + 0] = r;
		  		probeR[NprobeLocal*mesh.dim + 1] = s;
		  		probeR[NprobeLocal*mesh.dim + 2] = t;
		  		probeE[NprobeLocal] = e;
		  		probeIDl[NprobeLocal] = pid;
		   		probeB[p] = 1; // fix this prope
		   		NprobeLocal++; 
		   		break; 
				}else{
					const dfloat d1 = (r < -1.0) ? (-1.0 - r) : 0.0;
					const dfloat d2 = (s < -1.0) ? (-1.0 - s) : 0.0;
					const dfloat d3 = (t < -1.0) ? (-1.0 - t) : 0.0;
					const dfloat d4 = (r + s + t > -1.0) ? (r + s + t + 1.0) : 0.0;
					const dfloat dist2 =  0.25*(d1*d1 + d2*d2 + d3*d3 + d4*d4);

					const dfloat rc =bestR[p*mesh.dim + 0]; 
  				const dfloat sc =bestR[p*mesh.dim + 1]; 
  				const dfloat tc =bestR[p*mesh.dim + 2]; 

  				const dfloat dc1 = (rc < -1.0) ? (-1.0 - rc) : 0.0;
					const dfloat dc2 = (sc < -1.0) ? (-1.0 - sc) : 0.0;
					const dfloat dc3 = (tc < -1.0) ? (-1.0 - tc) : 0.0;
					const dfloat dc4 = (rc + sc + tc > -1.0) ? (rc + sc + tc + 1.0) : 0.0;
					const dfloat distc2 =  0.25*(dc1*dc1 + dc2*dc2 + dc3*dc3 + dc4*dc4);

					bestD[p] = dist2<distc2? dist2:distc2 ; 
		   	  bestR[p*mesh.dim + 0] = dist2<distc2? r: rc ; // r
			  	bestR[p*mesh.dim + 1] = dist2<distc2? s: sc ; // r
			  	bestR[p*mesh.dim + 2] = dist2<distc2? t: tc ; // r
			  	bestE[p]              = dist2<distc2? e: bestE[p];  

				}
				
		   }

		 }
		}
	}

void cns_t::reportProbe(const dfloat T, const dfloat tstep, int frame){

	std::string name;
	settings.getSetting("PROBE OUTPUT FILE", name);
	char fname[BUFSIZ];
	sprintf(fname, "%s.dat", name.c_str());  
	FILE *fp;

	if(mesh.rank==0){
		fp = fopen(fname, "a");
		if(frame==0){
			if(mesh.dim==2){
		fprintf(fp, "/* time probeID Pressure x-velocity y-Velocity */\n"); 
			}else{
		fprintf(fp, "/* time probeID Pressure x-velocity y-Velocity z-velocity*/\n");   		
			}  	
		}

	}

	const int eID = (mesh.dim==3) ? 4:3;

	memory<dfloat> p(mesh.Np*NprobeLocal);
	memory<dfloat> u(mesh.Np*NprobeLocal);
	memory<dfloat> v(mesh.Np*NprobeLocal);
	memory<dfloat> w(mesh.Np*NprobeLocal);

	memory<dfloat> Ip(NprobeLocal);
	memory<dfloat> Iu(NprobeLocal);
	memory<dfloat> Iv(NprobeLocal);
	memory<dfloat> Iw(NprobeLocal);


	for(int i=0; i<NprobeLocal; i++){
	// const int pid = probeID[p]; 
		const int e   = probeE[i]; 
		for(int n=0;n<mesh.Np;++n){
			dfloat rm = q[e*mesh.Np*Nfields+n];
			dfloat um = q[e*mesh.Np*Nfields+n+mesh.Np*1]/rm;
			dfloat vm = q[e*mesh.Np*Nfields+n+mesh.Np*2]/rm;
			dfloat wm = mesh.dim==3 ? q[e*mesh.Np*Nfields+n+mesh.Np*3]/rm:0.0;
			dfloat em = q[e*mesh.Np*Nfields+n+mesh.Np*eID];
			dfloat pm = (gamma-1)*(em-0.5*rm*(um*um+vm*vm+wm*wm));
	  //
			u[i*mesh.Np + n] = um; 
			v[i*mesh.Np + n] = vm; 
			w[i*mesh.Np + n] = wm; 
			p[i*mesh.Np + n] = pm; 

		}
	}

	probeInterp(p, Ip); 
	probeInterp(u, Iu); 
	probeInterp(v, Iv); 
	if(mesh.dim==3){ probeInterp(w, Iw);} 

	memory<dfloat> Ipg; 
	memory<dfloat> Iug; 
	memory<dfloat> Ivg; 
	memory<dfloat> Iwg; 

	if(mesh.rank==0){
		Ipg.malloc(NprobeGlobal);
		Iug.malloc(NprobeGlobal);
		Ivg.malloc(NprobeGlobal);
		Iwg.malloc(NprobeGlobal);
	}

	mesh.comm.Gatherv(Ip, NprobeLocal, Ipg, probeRecvCount,probeRecvOffset, 0);
	mesh.comm.Gatherv(Iu, NprobeLocal, Iug, probeRecvCount,probeRecvOffset, 0);
	mesh.comm.Gatherv(Iv, NprobeLocal, Ivg, probeRecvCount,probeRecvOffset, 0);
	if(mesh.dim==3){
		mesh.comm.Gatherv(Iw, NprobeLocal, Iwg, probeRecvCount, probeRecvOffset, 0);	
	}

	if(mesh.rank==0){
		for(int i=0; i<NprobeGlobal; i++){
			const int pid = probeIDa[i]; 
			if(mesh.dim==2){
				fprintf(fp, "%.6e %d %.6e %.6e %.6e\n", T, pid, Ipg[i],Iug[i], Ivg[i]);
			}else{ 
				fprintf(fp, "%.6e %d %.6e %.6e %.6e %.6e\n", T, pid, Ipg[i],Iug[i], Ivg[i],Iwg[i]);
			} 
		}
		fclose(fp);
	}


}


void cns_t::locateProbesQuad2D(){
	printf("Here quad locate function\n");
}

void cns_t::locateProbesHex3D(){
	printf("Here Hex locate function\n");
}








































// all positive barycentric coordinates
			// if(lmin>tol){
			// 	probeR[p*mesh.dim + 0] = 2.0*l3-1.0; // r
		  // 	probeR[p*mesh.dim + 1] = 2.0*l1-1.0; // s
		  //  	probeE[p] = e;
		  //  	probeB[p] = 1; // fix this prope
			// 	// inside = 1; 
			// 	// NprobeLocal++; 
			// 	break; 
			// }

			// // p is outside of the triangle
			// if(l1<0 || l2<0 || l3<0){
			// 	if(probeB[p]==0){
			// 		// printf("here for p = %d at e = %d \n", p, e);
			// 		dist2 = alpha*( (l1-alpha)*(l1-alpha)+(l2-alpha)*(l2-alpha)+(l3-alpha)*(l3-alpha)); 
			// 		if(dist2<bestdist2){
			// 			bestdist2 = dist2; 
			// 			probeR[p*mesh.dim + 0] = 2.0*l3-1.0; // r
			//   	  probeR[p*mesh.dim + 1] = 2.0*l1-1.0; // s
			//   	  probeE[p] = e;
			// 	   	probeB[p] = 0; // fix this prope
			// 			// NprobeLocal++; 
			// 		}
			// 	}
			// }



  	// Check whetver inside or outside







  	// if(e==18371){
  	// for(int p=0; p<NprobeGlobal; p++){
  	// 	printf("%d %d %.4e %.4e %.4e\n ",e, p, c[p*(mesh.dim+1) + 0], c[p*(mesh.dim+1) + 1], c[p*(mesh.dim+1) + 2]); 

  	// 	dfloat l1 = c[0*(mesh.dim+1) + 2]; 
  	// 	dfloat l2 = c[0*(mesh.dim+1) + 0]; 
  	// 	dfloat l3 = c[0*(mesh.dim+1) + 1]; 


  	// 	dfloat r = 2*l3-1.0; 
  	// 	dfloat s = 2*l1-1.0; 

  	// 	printf(" %.4e %.4e \n", r, s);



  	// }
  	// printf("-----------------------------\n");

  	// for(int p=0; p<NprobeGlobal; p++){




  // 	}



  // // }






	// }


// e= 6 18371 5.3737e-01 -9.0091e-01 
// e= 6 17490 -2.0244e-02 -2.1475e-01 
// e= 6 8600 -1.1591e-01 -3.9927e-04 
// e= 6 7698 -3.0006e-01 -2.9126e-01 
// e= 6 2466 -4.1735e-01 1.4906e-01 
// e= 6 1826 -8.4312e-01 -4.2532e-01 










// probeR.malloc(NprobeGlobal*mesh.dim); 
// probeE.malloc(NprobeGlobal); 
// NprobeLocal = 0; 

// for (int n = 0; n < NprobeGlobal; ++n) {
// 	dfloat bestdist2 = 1e12; 
// 	dfloat bestr     = 1e12; 
// 	dfloat bests     = 1e12; 
// 	int    beste     = -1; 
// 	int inside       = 0; 

// 	const dfloat px = probeX[n*mesh.dim+0];
// 	const dfloat py = probeX[n*mesh.dim+1];
// 	const int    pid= probeID[n]; 

// 	const dfloat tol = 1e-8;

// 	for(dlong e=0; e<mesh.Nelements; e++){
// 		const dfloat x1 = mesh.EX[e*mesh.Nverts+0];
// 		const dfloat x2 = mesh.EX[e*mesh.Nverts+1];
// 		const dfloat x3 = mesh.EX[e*mesh.Nverts+2];

// 		const dfloat y1 = mesh.EY[e*mesh.Nverts+0];
// 		const dfloat y2 = mesh.EY[e*mesh.Nverts+1];
// 		const dfloat y3 = mesh.EY[e*mesh.Nverts+2];
// 		// edge vectors
// 	  const dfloat xr = 0.5*(x2 - x1);
// 	  const dfloat xs = 0.5*(x3 - x1);

// 	  const dfloat yr = 0.5*(y2 - y1);
// 	  const dfloat ys = 0.5*(y3 - y1);

// 	  const dfloat pxs = px - x1 - xr - xs;
//     const dfloat pys = py - y1 - yr - ys;
// 		// --- Jacobian ---
// 		const dfloat J = xr*ys - xs*yr;
// 		const dfloat invJ = 1.0 / J;
// 		// --- solve for (r,s) ---
// 		dfloat r = invJ*( pxs*ys - pys*xs );
// 		dfloat s = invJ*(-pxs*yr + pys*xr );

// 		dfloat d1 = std::max(0.0, -1 - r);
// 		dfloat d2 = std::max(0.0, -1 - s);
// 		dfloat d3 = std::max(0.0,  r + s);

// 		dfloat dist2 = d1*d1 + d2*d2 + d3*d3;

// 		if(dist2<tol){
// 			bestdist2 = dist2; 
// 			bestr     = r; 
// 			bests     = s; 
// 			beste     = e; 
// 			inside    = 1; 
// 			break; 
// 		}
// 		if(dist2<bestdist2){
// 			bestdist2 = dist2; 
// 			bestr     = r; 
// 			bests     = s; 
// 			beste     = e; 
// 		}
// 	}

//       probeR[NprobeLocal*mesh.dim + 0] = bestr; 
// 			probeR[NprobeLocal*mesh.dim + 1] = bests;
// 			probeE[NprobeLocal] = beste;
// 			NprobeLocal++; 
//  }


//  probeR.realloc(NprobeLocal*mesh.dim); 
//  probeE.realloc(NprobeLocal); 

//  for(int n=0; n<NprobeLocal; n++){
//  	printf("e= %d %d %.4e %.4e \n", NprobeLocal, probeE[n], probeR[n*mesh.dim+0], probeR[n*mesh.dim+1]);
//  }
// }













// void cns_t::locateProbesTet3D(){

//   probeR.malloc(NprobeGlobal*mesh.dim); 
// 	probeE.malloc(NprobeGlobal); 

// 	NprobeLocal = 0; 

// 	for (int n = 0; n < NprobeGlobal; ++n) {
// 		const dfloat px = probeX[n*mesh.dim+0];
// 		const dfloat py = probeX[n*mesh.dim+1];
// 		const dfloat pz = probeX[n*mesh.dim+2];
// 		const int    pid= probeID[n]; 

// 		const dfloat tol = 1e-12;

// 		dfloat bestdist2 = 1e12; 
// 		dfloat bestr     = 1e12; 
// 		dfloat bests     = 1e12; 
// 		dfloat bestt     = 1e12; 
// 		int    beste     = -1; 

// 		int inside       = 0; 

// 		for(dlong e=0; e<mesh.Nelements; e++){
// 			const dfloat x1 = mesh.EX[e*mesh.Nverts+0];
// 			const dfloat x2 = mesh.EX[e*mesh.Nverts+1];
// 			const dfloat x3 = mesh.EX[e*mesh.Nverts+2];
// 			const dfloat x4 = mesh.EX[e*mesh.Nverts+3];

// 			const dfloat y1 = mesh.EY[e*mesh.Nverts+0];
// 			const dfloat y2 = mesh.EY[e*mesh.Nverts+1];
// 			const dfloat y3 = mesh.EY[e*mesh.Nverts+2];
// 			const dfloat y4 = mesh.EY[e*mesh.Nverts+3];

// 			const dfloat z1 = mesh.EZ[e*mesh.Nverts+0];
// 			const dfloat z2 = mesh.EZ[e*mesh.Nverts+1];
// 			const dfloat z3 = mesh.EZ[e*mesh.Nverts+2];
// 			const dfloat z4 = mesh.EZ[e*mesh.Nverts+3];	

// 					// edge vectors (HW geometric mapping)
// 		  const dfloat xr = 0.5*(x2 - x1);
// 		  const dfloat xs = 0.5*(x3 - x1);
// 		  const dfloat xt = 0.5*(x4 - x1);

// 		  const dfloat yr = 0.5*(y2 - y1);
// 		  const dfloat ys = 0.5*(y3 - y1);
// 		  const dfloat yt = 0.5*(y4 - y1);

// 		  const dfloat zr = 0.5*(z2 - z1);
// 		  const dfloat zs = 0.5*(z3 - z1);
// 		  const dfloat zt = 0.5*(z4 - z1);

// 		  // shift RHS 
// 		  const dfloat pxs = px - 0.5*(x2 + x3 + x4 - x1);
// 		  const dfloat pys = py - 0.5*(y2 + y3 + y4 - y1);
// 		  const dfloat pzs = pz - 0.5*(z2 + z3 + z4 - z1);

// 			  // Jacobian (constant)
// 			  const dfloat J =xr*(ys*zt - yt*zs)- xs*(yr*zt - yt*zr)+ xt*(yr*zs - ys*zr);


// 			  const dfloat invJ = 1.0 /(1.0*J);

// 			  // Solve using HW geometric factors
// 			  dfloat r = invJ*( pxs*(ys*zt - yt*zs)- pys*(xs*zt - xt*zs)+ pzs*(xs*yt - xt*ys) );
// 			  dfloat s = invJ*(-pxs*(yr*zt - yt*zr)+ pys*(xr*zt - xt*zr)- pzs*(xr*yt - xt*yr) );
// 			  dfloat t = invJ*( pxs*(yr*zs - ys*zr)- pys*(xr*zs - xs*zr)+ pzs*(xr*ys - xs*yr) );

// 			// bool inside = (r>=-1.0-tol)&&(s>=-1.0-tol)&&(t>=-1.0-tol) && (r+s+t<=-1.0+tol);
// 			// bool inside = (r>=-1.0)&&(s>=-1.0)&&(t>=-1.0) && (r+s+t<=-1.0);
// 			dfloat d1 = std::max(0.0, -1 - r);
// 			dfloat d2 = std::max(0.0, -1 - s);
// 			dfloat d3 = std::max(0.0, -1 - t);
// 			dfloat d4 = std::max(0.0,  r + s + t + 1);

// 			dfloat dist2 = d1*d1 + d2*d2 + d3*d3 + d4*d4;

// 			if(dist2<tol){
// 				bestdist2 = dist2; 
// 				bestr     = r; 
// 				bests     = s; 
// 				bestt     = t; 
// 				beste     = e; 
// 				inside    = 1; 
// 				break; 
// 			}


// 			if(dist2<bestdist2){
// 				bestdist2 = dist2; 
// 				bestr     = r; 
// 				bests     = s; 
// 				bestt     = t; 
// 				beste     = e; 
// 			}

// 			// if(inside){
// 			// 	probeR[NprobeLocal*mesh.dim + 0] = r; 
// 			// 	probeR[NprobeLocal*mesh.dim + 1] = s;
// 			// 	probeR[NprobeLocal*mesh.dim + 2] = t;
// 			// 	probeE[NprobeLocal] = e;
// 			// 	NprobeLocal++; 
// 			// }
// 		}

//         probeR[NprobeLocal*mesh.dim + 0] = bestr; 
// 				probeR[NprobeLocal*mesh.dim + 1] = bests;
// 				probeR[NprobeLocal*mesh.dim + 2] = bestt;
// 				probeE[NprobeLocal] = beste;
// 				NprobeLocal++; 

// 				printf("beste = %d %d %.4e, %d\n", beste, inside, bestdist2, pid); 
// 	 }

// 	 probeR.realloc(NprobeLocal*mesh.dim); 
// 	 probeE.realloc(NprobeLocal); 

// 	 for(int n=0; n<NprobeLocal; n++){
// 	 	printf("e= %d %d %.4e %.4e %.4e\n", 
// 	 		      NprobeLocal, probeE[n], probeR[n*mesh.dim+0], probeR[n*mesh.dim+1],  probeR[n*mesh.dim+2]);
// 	 }


// }


