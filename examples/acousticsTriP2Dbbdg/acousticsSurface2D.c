#include "acoustics2D.h"

void boundaryConditions2D(iint bc, dfloat t, dfloat x, dfloat y,
			dfloat uM, dfloat vM, dfloat pM,
						dfloat *uP, dfloat *vP, dfloat *pP){
	if(1){//bc==1){
		*uP = -uM;  
		*vP = -vM;  
		*pP = pM; 
	}   
	if(0){ // (bc==2){  
		dfloat dx = 1.f/sqrt(2.f);
		dfloat dy = 1.f/sqrt(2.f);
		dfloat omega = 10.f*M_PI;
		dfloat wave = cos(omega*(t-(x*dx+y*dy))); 
		dfloat uI = dx*wave;
		dfloat vI = dy*wave;
		dfloat pI = wave; 
		*uP = -uM -2.f*uI;  
		*vP = -vM -2.f*vI;
		*pP = pM;   
	}
}

// function to compute surface contributions 
// for nodal DG acoustics right hand side
void acousticsSurface2Dbbdg(mesh2D *mesh, dfloat t){

	// temporary storage for flux terms
	dfloat *fluxu = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
	dfloat *fluxv = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
	dfloat *fluxp = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

	dfloat *fluxu_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
	dfloat *fluxv_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));
	dfloat *fluxp_copy = (dfloat*) calloc(mesh->NfpMax*mesh->Nfaces,sizeof(dfloat));

	dfloat nx, ny, sJ, invJ;

	// for all elements
	for(iint e=0;e<mesh->Nelements;++e){
		iint N = mesh->N[e];
		for (iint f=0;f<mesh->Nfaces;f++){
			// load element and face number of neighbour
			iint eP = mesh->EToE[e*mesh->Nfaces+f];
			iint fP = mesh->EToF[e*mesh->Nfaces+f];
			iint NP = mesh->N[eP]; 
			
			if (eP<0 || fP<0) NP = N; //boundary

			for (iint n=0;n<mesh->Nfp[NP];n++){
				iint id  = e*mesh->Nfaces*mesh->NfpMax + f*mesh->NfpMax + n;
				iint idP = mesh->vmapP[id];
				iint qidP = mesh->Nfields*idP;

				//load qP into flux for now to save space
				fluxu[n] = mesh->q[qidP+0];
				fluxv[n] = mesh->q[qidP+1];
				fluxp[n] = mesh->q[qidP+2];
			}

			if (NP < N) { 
				for (iint n=0;n<mesh->Nfp[N];n++){
					fluxu_copy[f*mesh->Nfp[N] + n] = 0.0;
					fluxv_copy[f*mesh->Nfp[N] + n] = 0.0;
					fluxp_copy[f*mesh->Nfp[N] + n] = 0.0;
					for (iint m=0;m<2;m++){ //apply raise operator sparsly
						dfloat BBRaiseVal = mesh->BBRaiseVals[N][2*n+m];
						iint BBRaiseid = mesh->BBRaiseids[N][2*n+m];
						fluxu_copy[f*mesh->Nfp[N] + n] += BBRaiseVal*fluxu[BBRaiseid];
						fluxv_copy[f*mesh->Nfp[N] + n] += BBRaiseVal*fluxv[BBRaiseid];
						fluxp_copy[f*mesh->Nfp[N] + n] += BBRaiseVal*fluxp[BBRaiseid];
					}
				}
			} else if (NP > N) { 
				for (iint n=0;n<mesh->Nfp[N];n++){
					fluxu_copy[f*mesh->Nfp[N] +n] = 0.0;
					fluxv_copy[f*mesh->Nfp[N] +n] = 0.0;
					fluxp_copy[f*mesh->Nfp[N] +n] = 0.0;
					for (iint m=0;m<mesh->Nfp[NP];m++){
						iint id = n*mesh->Nfp[NP] + m;
						fluxu_copy[f*mesh->Nfp[N] + n] += mesh->BBLower[N][id]*fluxu[m];
						fluxv_copy[f*mesh->Nfp[N] + n] += mesh->BBLower[N][id]*fluxv[m];
						fluxp_copy[f*mesh->Nfp[N] + n] += mesh->BBLower[N][id]*fluxp[m];
					}
				}
			} else { //equal order neighbor
				for (iint n=0;n<mesh->Nfp[N];n++){
					//load qP into flux_copy to save space
					fluxu_copy[f*mesh->Nfp[N] + n] = fluxu[n];
					fluxv_copy[f*mesh->Nfp[N] + n] = fluxv[n];
					fluxp_copy[f*mesh->Nfp[N] + n] = fluxp[n];
				}
			}
		}
		

		// for all face nodes of all elements
		for (iint f=0;f<mesh->Nfaces;f++) {
			for(iint n=0;n<mesh->Nfp[N];++n){
				iint eP = mesh->EToE[e*mesh->Nfaces+f];
				iint fP = mesh->EToF[e*mesh->Nfaces+f];
				iint NP = mesh->N[eP];
				if (eP <0 || fP<0) NP = N;

				// load surface geofactors for this face
				iint sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
				nx   = mesh->sgeo[sid+0];
				ny   = mesh->sgeo[sid+1];
				sJ   = mesh->sgeo[sid+2];
				invJ = mesh->sgeo[sid+3];

				// indices of negative and positive traces of face node
				iint id  = e*mesh->NfpMax*mesh->Nfaces +f*mesh->NfpMax + n;
				iint idM = mesh->vmapM[id];
				iint idP = mesh->vmapP[id];
				iint qidM = mesh->Nfields*idM;

				// load negative trace node values of q
				dfloat uM = mesh->q[qidM+0];
				dfloat vM = mesh->q[qidM+1];
				dfloat pM = mesh->q[qidM+2];
				
				// load positive trace node values of q from flux_copy
				//dfloat uP = mesh->q[qidP+0];
				//dfloat vP = mesh->q[qidP+1];
				//dfloat pP = mesh->q[qidP+2];
				dfloat uP = fluxu_copy[f*mesh->Nfp[N] + n]; 
				dfloat vP = fluxv_copy[f*mesh->Nfp[N] + n];
				dfloat pP = fluxp_copy[f*mesh->Nfp[N] + n];

				// find boundary type
				iint boundaryType = mesh->EToB[e*mesh->Nfaces+f];
				if(boundaryType>0)
				boundaryConditions2D(boundaryType, t, mesh->x[idM], mesh->y[idM], uM, vM, pM, &uP, &vP, &pP);

				// compute (q^* - q^-)
				dfloat duS = 0.5f*(uP-uM) + mesh->Lambda2*(-nx)*(pP-pM);
				dfloat dvS = 0.5f*(vP-vM) + mesh->Lambda2*(-ny)*(pP-pM);
				dfloat dpS = 0.5f*(pP-pM) + mesh->Lambda2*(-nx*(uP-uM)-ny*(vP-vM));

				// evaluate "flux" terms: (sJ/J)*(A*nx+B*ny)*(q^* - q^-)
				fluxu[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*dpS);
				fluxv[f*mesh->Nfp[N] +n] = invJ*sJ*(-ny*dpS);
				fluxp[f*mesh->Nfp[N] +n] = invJ*sJ*(-nx*duS-ny*dvS);
			}
		}

		// apply L0 to fluxes. use fact that L0 = tridiagonal in 2D
		for(iint n=0;n<mesh->Nfp[N]*mesh->Nfaces;++n){
		
			iint id = n % mesh->Nfp[N];  // warning: redundant reads
			dfloat L0val = mesh->L0vals[N][3*id+1]; 

			dfloat utmpflux = L0val * fluxu[n];
			dfloat vtmpflux = L0val * fluxv[n];
			dfloat ptmpflux = L0val * fluxp[n];

			if (id > 0){     
				utmpflux += mesh->L0vals[N][3*id]*fluxu[n-1]; // add previous term
				vtmpflux += mesh->L0vals[N][3*id]*fluxv[n-1]; // add previous term
				ptmpflux += mesh->L0vals[N][3*id]*fluxp[n-1]; // add previous term
			}
			if (id < mesh->Nfp[N]){
				utmpflux += mesh->L0vals[N][3*id+2]*fluxu[n+1];// add next term
				vtmpflux += mesh->L0vals[N][3*id+2]*fluxv[n+1];// add next term
				ptmpflux += mesh->L0vals[N][3*id+2]*fluxp[n+1];// add next term
			}
			fluxu_copy[n] = utmpflux;
			fluxv_copy[n] = vtmpflux;
			fluxp_copy[n] = ptmpflux;
		}

		// apply lift reduction and accumulate RHS
    for(iint n=0;n<mesh->Np[N];++n){
	    iint id = mesh->Nfields*(mesh->NpMax*e + n);
	    
	    // load RHS
	    dfloat rhsqnu = mesh->rhsq[id+0];
	    dfloat rhsqnv = mesh->rhsq[id+1];
	    dfloat rhsqnp = mesh->rhsq[id+2];

	    // rhs += LIFT*((sJ/J)*(A*nx+B*ny)*(q^* - q^-))
      for (int m = 0; m < mesh->max_EL_nnz[N]; ++m){
				iint id = m + n*mesh->max_EL_nnz[N];
				dfloat ELval = mesh->ELvals[N][id];
				iint ELid = mesh->ELids[N][id];
				rhsqnu += ELval * fluxu_copy[ELid];
				rhsqnv += ELval * fluxv_copy[ELid];
				rhsqnp += ELval * fluxp_copy[ELid];
      }
	    
	    // store incremented rhs
	    mesh->rhsq[id+0] = rhsqnu;
	    mesh->rhsq[id+1] = rhsqnv;
	    mesh->rhsq[id+2] = rhsqnp;  
    }
  }

	// free temporary flux storage
	free(fluxu); free(fluxv); free(fluxp);
	free(fluxu_copy); free(fluxv_copy); free(fluxp_copy);
}
