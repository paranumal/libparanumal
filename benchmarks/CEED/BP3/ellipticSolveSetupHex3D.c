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

#include "ellipticHex3D.h"

occa::kernel saferBuildKernelFromSource(occa::device &device,
                                        const char *fname, const char *kname, occa::kernelInfo &kernelInfo){
                                        
	// should really use root to build and non-root to load
	return device.buildKernelFromSource(fname, kname, kernelInfo);
	
}



// specialized version for geometric factors at Gauss (not GLL) nodes
//dfloat *ellipticGeometricFactorsHex3D(mesh3D *mesh){
dfloat  *ellipticGeometricFactorsHex3D(mesh3D *mesh){
	/* number of second order geometric factors */
	int NgjGeo = 7;
	int gjNq = mesh->gjNq;
	int gjNp = gjNq*gjNq*gjNq;
	dfloat *gjGeo = (dfloat*) calloc(mesh->Nelements*NgjGeo*gjNp, sizeof(dfloat));
	//KS start
	printf("Setup info: Nggeo = %d, gjNq = %d. gjNp = %d Nq = %d \n", NgjGeo, gjNq, gjNp, mesh->Nq);
	#if 1
	//dfloat *gjGeod2 = (dfloat*) calloc(mesh->Nelements*8*gjNp, sizeof(dfloat));
	#else
	//dfloat *gjGeod2 = (dfloat*) calloc(mesh->Nelements*3*8+2*gjNp, sizeof(dfloat));
	#endif
	//KS end
	for(int e=0; e<mesh->Nelements; ++e) { /* for each element */
	
		/* find vertex indices and physical coordinates */
		int id = e*mesh->Nverts;
		
		dfloat *xe = mesh->EX + id;
		dfloat *ye = mesh->EY + id;
		dfloat *ze = mesh->EZ + id;
		
		for(int k=0; k<gjNq; ++k) {
			for(int j=0; j<gjNq; ++j) {
				for(int i=0; i<gjNq; ++i) {
				
					int n = i + j*gjNq + k*gjNq*gjNq;
					
					/* local node coordinates */
					dfloat rn = mesh->gjr[i];
					dfloat sn = mesh->gjr[j];
					dfloat tn = mesh->gjr[k];
					//	printf("r,s,t=%g,%g,%g\n", rn,sn,tn);
					/* Jacobian matrix */
					dfloat xr = 0.125*( (1-tn)*(1-sn)*(xe[1]-xe[0]) + (1-tn)*(1+sn)*(xe[2]-xe[3]) + (1+tn)*(1-sn)*(xe[5]-xe[4]) + (1+tn)*(1+sn)*(xe[6]-xe[7]) );
					dfloat xs = 0.125*( (1-tn)*(1-rn)*(xe[3]-xe[0]) + (1-tn)*(1+rn)*(xe[2]-xe[1]) + (1+tn)*(1-rn)*(xe[7]-xe[4]) + (1+tn)*(1+rn)*(xe[6]-xe[5]) );
					dfloat xt = 0.125*( (1-rn)*(1-sn)*(xe[4]-xe[0]) + (1+rn)*(1-sn)*(xe[5]-xe[1]) + (1+rn)*(1+sn)*(xe[6]-xe[2]) + (1-rn)*(1+sn)*(xe[7]-xe[3]) );
					
					dfloat yr = 0.125*( (1-tn)*(1-sn)*(ye[1]-ye[0]) + (1-tn)*(1+sn)*(ye[2]-ye[3]) + (1+tn)*(1-sn)*(ye[5]-ye[4]) + (1+tn)*(1+sn)*(ye[6]-ye[7]) );
					dfloat ys = 0.125*( (1-tn)*(1-rn)*(ye[3]-ye[0]) + (1-tn)*(1+rn)*(ye[2]-ye[1]) + (1+tn)*(1-rn)*(ye[7]-ye[4]) + (1+tn)*(1+rn)*(ye[6]-ye[5]) );
					dfloat yt = 0.125*( (1-rn)*(1-sn)*(ye[4]-ye[0]) + (1+rn)*(1-sn)*(ye[5]-ye[1]) + (1+rn)*(1+sn)*(ye[6]-ye[2]) + (1-rn)*(1+sn)*(ye[7]-ye[3]) );
					
					dfloat zr = 0.125*( (1-tn)*(1-sn)*(ze[1]-ze[0]) + (1-tn)*(1+sn)*(ze[2]-ze[3]) + (1+tn)*(1-sn)*(ze[5]-ze[4]) + (1+tn)*(1+sn)*(ze[6]-ze[7]) );
					dfloat zs = 0.125*( (1-tn)*(1-rn)*(ze[3]-ze[0]) + (1-tn)*(1+rn)*(ze[2]-ze[1]) + (1+tn)*(1-rn)*(ze[7]-ze[4]) + (1+tn)*(1+rn)*(ze[6]-ze[5]) );
					dfloat zt = 0.125*( (1-rn)*(1-sn)*(ze[4]-ze[0]) + (1+rn)*(1-sn)*(ze[5]-ze[1]) + (1+rn)*(1+sn)*(ze[6]-ze[2]) + (1-rn)*(1+sn)*(ze[7]-ze[3]) );
					
					/* compute geometric factors for affine coordinate transform*/
					dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
					
					if(J<1e-12) printf("J = %g !!!!!!!!!!!!!\n", J);
					
					dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
					dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
					dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;
					
					dfloat JW = J*mesh->gjw[i]*mesh->gjw[j]*mesh->gjw[k];
					
					/* store second order geometric factors */
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
					gjGeo[NgjGeo*gjNp*e + n + gjNp*GWJID] = JW;
					
					
					//emap*p_Nggeo*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;
					//ggeo[base+p_Np*4+0];
					/*			gjGeod2[8*gjNp*e + 8*n + 0 ] = JW*(rx*rx + ry*ry + rz*rz);
								gjGeod2[8*gjNp*e + 8*n + 1 ] = JW*(rx*sx + ry*sy + rz*sz);
								
								gjGeod2[8*gjNp*e + 8*n + 2 ] = JW*(rx*tx + ry*ty + rz*tz);
								gjGeod2[8*gjNp*e + 8*n + 3 ] = JW*(sx*sx + sy*sy + sz*sz);
								
								gjGeod2[8*gjNp*e + 8*n + 4 ] = JW*(sx*tx + sy*ty + sz*tz);
								gjGeod2[8*gjNp*e + 8*n + 5 ] = JW*(tx*tx + ty*ty + tz*tz);
								
								gjGeod2[8*gjNp*e + 8*n + 6 ] = JW;
								gjGeod2[8*gjNp*e + 8*n + 7 ] = 0.0;
								*/
					/*	  gjGeod2[8*gjNp*e + n + gjNp*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
							gjGeod2[8*gjNp*e + n + gjNp*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
							gjGeod2[8*gjNp*e + n + gjNp*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
							gjGeod2[8*gjNp*e + n + gjNp*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
							gjGeod2[8*gjNp*e + n + gjNp*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
							gjGeod2[8*gjNp*e + n + gjNp*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
							gjGeod2[8*gjNp*e + n + gjNp*GWJID] = JW;
							gjGeod2[8*gjNp*e + n + gjNp*7] = 0.0f;*/
					
					
				}
			}
		}
	}
	/*
	#if 1
	for (int n = 0; n< mesh->Nelements*8*gjNp; n++)
{
		gjGeod2[n] = drand48();;
}
	#else
	for(int k=0; k<gjNq; k++) {
			// per element, load coordinated and weights
			gjGeod2[2*k ]  = mesh->gjr[k];
			gjGeod2[2*k +1]  = mesh->gjw[k];
		}
	
	for(int e=0; e<mesh->Nelements; ++e) {
		
		int id = e*mesh->Nverts;
		for (int k= 0; k<8; k++)
		{
			// in the memory: first all x coordinates for all the elements
			// then y coordinates
			// then z coordinates 
	
			gjGeod2[2*gjNq + e*mesh->Nelements + k] = mesh->EX[id+k]; 					
			gjGeod2[2*gjNq + mesh->Nelements*8 + e*mesh->Nelements + k] = mesh->EY[id+k];		
			gjGeod2[2*gjNq + mesh->Nelements*16+ e*mesh->Nelements + k] = mesh->EZ[id+k];
	
		} 
	
	
}

	#endif
	*/
	
	return gjGeo;
}


void ellipticComputeDegreeVector(mesh3D *mesh, int Ntotal, ogs_t *ogs, dfloat *deg){

	// build degree vector
	for(int n=0; n<Ntotal; ++n)
		deg[n] = 1;
		
	occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
	
	o_deg.copyFrom(deg);
	
	ellipticParallelGatherScatter(mesh, ogs, o_deg, o_deg, dfloatString, "add");
	
	o_deg.copyTo(deg);
	
	mesh->device.finish();
	o_deg.free();
	
}

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options) {

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
	int NblockV = mymax(1,1024/mesh->Np); // works for CUDA
	int NblockV2 = 1; // works for CUDA
	int NblockS = mymax(1,1024/maxNodes); // works for CUDA
	int NblockG;
	
	int gjNq = mesh->gjNq;
	int gjNp = gjNq*gjNq*gjNq;
	int gjNq2 = gjNq*gjNq;
	if(gjNq2<=32)
		NblockG = ( 32/gjNq2 );
	else {
		if(mesh->Nq<=6) {
			NblockG = 256/gjNq2;
		}
		else
			NblockG = 1;
	}
	//  NblockG = 512/gNq2;
	
	int Ntotal = mesh->Np*mesh->Nelements;
	int NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;
	printf("NqP = %d \n", mesh->NqP);
	int Nblock = (Ntotal+blockSize-1)/blockSize;
	
	int Nhalo = mesh->Np*mesh->totalHaloPairs;
	int Nall   = Ntotal + Nhalo;
	int NallP  = NtotalP;
	
	solver_t *solver = (solver_t*) calloc(1, sizeof(solver_t));
	
	solver->mesh = mesh;
	
	solver->p   = (dfloat*) calloc(Nall, sizeof(dfloat));
	solver->r   = (dfloat*) calloc(Nall, sizeof(dfloat));
	solver->z   = (dfloat*) calloc(Nall, sizeof(dfloat));
	solver->zP  = (dfloat*) calloc(NallP, sizeof(dfloat));
	solver->Ax  = (dfloat*) calloc((mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nq+1)*mesh->Nelements, sizeof(dfloat));
	
	
	solver->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
	//	solver->grad = (dfloat*) calloc(4*(Ntotal+Nhalo), sizeof(dfloat));
	solver->grad = (dfloat*) calloc(4*((mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nq+1)*mesh->Nelements+Nhalo), sizeof(dfloat));
	solver->o_p   = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_rtmp= mesh->device.malloc(Nall*sizeof(dfloat), solver->r);
	solver->o_z   = mesh->device.malloc(Nall*sizeof(dfloat), solver->z);
	solver->o_zP  = mesh->device.malloc(NallP*sizeof(dfloat), solver->zP);// CAUTION
	solver->o_Ax  = mesh->device.malloc((mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nq+1)*mesh->Nelements*sizeof(dfloat), solver->Ax);
	solver->o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), solver->Ap);
	solver->o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), solver->tmp);
	solver->o_grad  = mesh->device.malloc(4*((mesh->Nq+1)*(mesh->Nq+1)*(mesh->Nq+1)*mesh->Nelements +Nhalo )*sizeof(dfloat), solver->grad);
	solver->o_pAp  = mesh->device.malloc(sizeof(dfloat));
	
	solver->o_Aw  = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_w    = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	solver->o_s    = mesh->device.malloc(Nall*sizeof(dfloat), solver->p);
	
	int Nbytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
	
#if 0
	solver->sendBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
	solver->recvBuffer = (dfloat*) calloc(Nbytes/sizeof(dfloat), sizeof(dfloat));
#else
	solver->defaultStream = mesh->device.getStream();
	solver->dataStream = mesh->device.createStream();
	mesh->device.setStream(solver->defaultStream);
	
	if(Nbytes>0) {
		occa::memory o_sendBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
		occa::memory o_recvBuffer = mesh->device.mappedAlloc(Nbytes, NULL);
		
		solver->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
		solver->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
	}else{
		solver->sendBuffer = NULL;
		solver->recvBuffer = NULL;
	}
	mesh->device.setStream(solver->defaultStream);
#endif
	solver->Nblock = Nblock;
	
	// BP3 specific stuff starts here
	
	//  dfloat *gjGeo = ellipticGeometricFactorsHex3D(mesh);
	dfloat *gjGeo = ellipticGeometricFactorsHex3D(mesh);
#if 0
	// build a gjD that maps gjNq to gjNq
	dfloat *gjD2 = (dfloat*) calloc(gjNq*gjNq, sizeof(dfloat));
	for(int n=0;n<gjNq;++n){
		for(int m=0;m<gjNq;++m){
			for(int i=0;i<mesh->Nq;++i){
				gjD2[n*gjNq+m] += mesh->gjD[n*mesh->Nq+i]*mesh->gjI[m*mesh->Nq+i];
				
			}
		}
	}
#endif
	
	// TW: temporarily resize gjD
	mesh->gjD = (dfloat*) realloc(mesh->gjD, gjNq*gjNq*sizeof(dfloat));
	for (int i=0; i<gjNq; ++i){
		for (int j=0; j<gjNq; ++j){
			printf(" %16.16f ",i,j, mesh->gjD[i*gjNq+j]);
		}
		printf("\n");
	}
	solver->o_gjD = mesh->device.malloc(gjNq*gjNq*sizeof(dfloat), mesh->gjD);
	solver->o_gjD2 = mesh->device.malloc(gjNq*gjNq*sizeof(dfloat), mesh->gjD2);
	solver->o_gjI = mesh->device.malloc(gjNq*mesh->Nq*sizeof(dfloat), mesh->gjI);
	// solver->o_gjGeo = mesh->device.malloc(mesh->Nggeo*gjNp*mesh->Nelements*sizeof(dfloat), gjGeo);
	//KS
	solver->o_gjGeo = mesh->device.malloc(8*gjNp*mesh->Nelements*sizeof(dfloat), gjGeo);
	// BP3 specific stuff ends here
	
	kernelInfo.addParserFlag("automate-add-barriers", "disabled");
	
	kernelInfo.addCompilerFlag("-Xptxas -dlcm=ca");
	//  kernelInfo.addCompilerFlag("-G");
	kernelInfo.addCompilerFlag("-O3");
	
	// generically used for blocked DEVICE reductions
	kernelInfo.addDefine("p_blockSize", blockSize);
	
	printf("p_blockSize = %d \n", blockSize);
	kernelInfo.addDefine("p_maxNodes", maxNodes);
	kernelInfo.addDefine("p_Nmax", maxNodes);
	
	//kernelInfo.addDefine("p_halfD", (int) (mesh->Nq+mesh->Nq%2)/2);
	kernelInfo.addDefine("p_halfD", (mesh->gjNq+1)/2);
	printf("p_halfD = %f \n", (float) (mesh->gjNq+1)/2.0);
	kernelInfo.addDefine("p_NblockV", NblockV);
	kernelInfo.addDefine("p_NblockV2", NblockV2);
	kernelInfo.addDefine("p_NblockS", NblockS);
	
	kernelInfo.addDefine("p_NblockG", NblockG);
	
	kernelInfo.addDefine("p_Lambda2", 0.5f);
	
	kernelInfo.addDefine("p_gjNq", mesh->gjNq);
	kernelInfo.addDefine("p_NqP", (mesh->Nq+2));
	kernelInfo.addDefine("p_NpP", (mesh->NqP*mesh->NqP*mesh->NqP));
	kernelInfo.addDefine("p_Nverts", mesh->Nverts);
	
	int Nz = mymin(mesh->Nq, 64/mesh->Nq);
	kernelInfo.addDefine("p_Nz", Nz);
	printf("Nz = %d\n", Nz);
	
	//  occa::setVerboseCompilation(0);
	
	for(int r=0;r<size;++r){
		MPI_Barrier(MPI_COMM_WORLD);
		if(r==rank){
			printf("Building kernels for rank %d\n", rank);
			fflush(stdout);
			mesh->haloExtractKernel =
			  saferBuildKernelFromSource(mesh->device,
			                             LIBP_DIR "/okl/meshHaloExtract3D.okl",
			                             "meshHaloExtract3D",
			                             kernelInfo);
			                             
			mesh->gatherKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/gather.okl",
			                             "gather",
			                             kernelInfo);
			                             
			mesh->scatterKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/scatter.okl",
			                             "scatter",
			                             kernelInfo);
			                             
			mesh->gatherScatterKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/gatherScatter.okl",
			                             "gatherScatter",
			                             kernelInfo);
			                             
			                             
			mesh->getKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/get.okl",
			                             "get",
			                             kernelInfo);
			                             
			mesh->putKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/put.okl",
			                             "put",
			                             kernelInfo);
			                             
#if 0
			// WARNING
			if(mesh->Nq<12) {
				solver->AxKernel =
				  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticAxHex3D.okl",
				                             "ellipticAxHex3D_e3",
				                             kernelInfo);
			}
#endif
			
			// CPU version is e8, GPU version is e6
			// KS
			
			printf("building e9 kernel, halfD = %d \n", (mesh->Nq+mesh->Nq%2)/2);
			solver->partialAxKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticAxHex3DNOCOLLOCATIONEx.okl",
			                             "ellipticAxHex3D_Ref10",
			                             //"ellipticAxHex3D_cuboid0",
			                             //						     "ellipticAxHex3D_cube1",
			                             kernelInfo);
			                             
			                             
			mesh->weightedInnerProduct1Kernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/weightedInnerProduct1.okl",
			                             "weightedInnerProduct1",
			                             kernelInfo);
			                             
			saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/weightedInnerProduct2.okl",
			                           "weightedInnerProduct2",
			                           kernelInfo);
			                           
			mesh->innerProductKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/innerProduct.okl",
			                             "innerProduct",
			                             kernelInfo);
			                             
			mesh->scaledAddKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/scaledAdd.okl",
			                             "scaledAdd",
			                             kernelInfo);
			                             
			mesh->dotMultiplyKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/dotMultiply.okl",
			                             "dotMultiply",
			                             kernelInfo);
			                             
			mesh->dotDivideKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/dotDivide.okl",
			                             "dotDivide",
			                             kernelInfo);
			                             
			solver->gradientKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticGradientHex3D.okl",
			                             "ellipticGradientHex3D",
			                             kernelInfo);
			                             
			solver->partialGradientKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticGradientHex3D.okl",
			                             "ellipticPartialGradientHex3D",
			                             kernelInfo);
			                             
			if(mesh->Nq<12) {
				solver->ipdgKernel =
				  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticAxIpdgHex3D.okl",
				                             "ellipticAxIpdgHex3D",
				                             kernelInfo);
				                             
				solver->partialIpdgKernel =
				  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticAxIpdgHex3D.okl",
				                             "ellipticPartialAxIpdgHex3D",
				                             kernelInfo);
			}
			
			solver->combinedInnerProductKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticCombinedInnerProduct.okl",
			                             "ellipticCombinedInnerProduct",
			                             kernelInfo);
			                             
			solver->combinedUpdateKernel =
			  saferBuildKernelFromSource(mesh->device, LIBP_DIR "/okl/ellipticCombinedUpdate.okl",
			                             "ellipticCombinedUpdate",
			                             kernelInfo);
			usleep(8000);
		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	occaTimerTic(mesh->device,"GatherScatterSetup");
	
	// set up gslib MPI gather-scatter and OCCA gather/scatter arrays
	solver->ogs = meshParallelGatherScatterSetup(mesh,
	              mesh->Np*mesh->Nelements,
	              sizeof(dfloat),
	              mesh->gatherLocalIds,
	              mesh->gatherBaseIds,
	              mesh->gatherHaloFlags);
	occaTimerToc(mesh->device,"GatherScatterSetup");
	
	occaTimerTic(mesh->device,"DegreeVectorSetup");
	dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
	dfloat *degree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
	
	solver->o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);
	
	ellipticComputeDegreeVector(mesh, Ntotal, solver->ogs, degree);
	
	for(int n=0; n<Ntotal; ++n) { // need to weight inner products{
		if(degree[n] == 0) printf("WARNING!!!!\n");
		invDegree[n] = 1./degree[n];
	}
	
	solver->o_invDegree.copyFrom(invDegree);
	occaTimerToc(mesh->device,"DegreeVectorSetup");
	
	//fill geometric factors in halo
	if(mesh->totalHaloPairs) {
		int Nlocal = mesh->Nelements*mesh->Np;
		int Nhalo = mesh->totalHaloPairs*mesh->Np;
		
		dfloat *vgeoSendBuffer = (dfloat*) calloc(Nhalo*mesh->Nvgeo, sizeof(dfloat));
		
		// import geometric factors from halo elements
		mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat));
		
		meshHaloExchange(mesh,
		                 mesh->Nvgeo*mesh->Np*sizeof(dfloat),
		                 mesh->vgeo,
		                 vgeoSendBuffer,
		                 mesh->vgeo + Nlocal*mesh->Nvgeo);
		                 
		mesh->o_vgeo =
		  mesh->device.malloc((Nlocal+Nhalo)*mesh->Nvgeo*sizeof(dfloat), mesh->vgeo);
	}
	
	// build weights for continuous SEM L2 project --->
	dfloat *localMM = (dfloat*) calloc(Ntotal, sizeof(dfloat));
	
	for(int e=0; e<mesh->Nelements; ++e) {
		for(int n=0; n<mesh->Np; ++n) {
			dfloat wJ = mesh->ggeo[e*mesh->Np*mesh->Nggeo + n + GWJID*mesh->Np];
			localMM[n+e*mesh->Np] = wJ;
		}
	}
	
	occa::memory o_localMM = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
	occa::memory o_MM      = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
	
	// sum up all contributions at base nodes and scatter back
	
	ellipticParallelGatherScatter(mesh, solver->ogs, o_localMM, o_MM, dfloatString, "add");
	
	mesh->o_projectL2 = mesh->device.malloc(Ntotal*sizeof(dfloat), localMM);
	mesh->dotDivideKernel(Ntotal, o_localMM, o_MM, mesh->o_projectL2);
	
	free(localMM); o_MM.free(); o_localMM.free();
	
	if(rank==0)
		printf("starting elliptic parallel gather scatter setup\n");
		
	// set up separate gather scatter infrastructure for halo and non halo nodes
	//  mesh->device.setStream(solver->dataStream);
	ellipticParallelGatherScatterSetup(mesh,
	                                   mesh->Np*mesh->Nelements,
	                                   sizeof(dfloat),
	                                   mesh->gatherLocalIds,
	                                   mesh->gatherBaseIds,
	                                   mesh->gatherHaloFlags,
	                                   &(solver->halo),
	                                   &(solver->nonHalo));
	//  mesh->device.setStream(solver->defaultStream);
	
	
	// count elements that contribute to global C0 gather-scatter
	int globalCount = 0;
	int localCount = 0;
	int *localHaloFlags = (int*) calloc(mesh->Np*mesh->Nelements, sizeof(int));
	
	for(int n=0; n<mesh->Np*mesh->Nelements; ++n) {
		localHaloFlags[mesh->gatherLocalIds[n]] += mesh->gatherHaloFlags[n];
	}
	
	for(int e=0; e<mesh->Nelements; ++e) {
		int isHalo = 0;
		for(int n=0; n<mesh->Np; ++n) {
			if(localHaloFlags[e*mesh->Np+n]>0) {
				isHalo = 1;
			}
			if(localHaloFlags[e*mesh->Np+n]<0) {
				printf("found halo flag %d\n", localHaloFlags[e*mesh->Np+n]);
			}
		}
		globalCount += isHalo;
		localCount += 1-isHalo;
	}
	
	//  printf("local = %d, global = %d\n", localCount, globalCount);
	
	solver->globalGatherElementList    = (int*) calloc(globalCount, sizeof(int));
	solver->localGatherElementList = (int*) calloc(localCount, sizeof(int));
	
	globalCount = 0;
	localCount = 0;
	
	for(int e=0; e<mesh->Nelements; ++e) {
		int isHalo = 0;
		for(int n=0; n<mesh->Np; ++n) {
			if(localHaloFlags[e*mesh->Np+n]>0) {
				isHalo = 1;
			}
		}
		if(isHalo) {
			solver->globalGatherElementList[globalCount++] = e;
		}
		else{
			solver->localGatherElementList[localCount++] = e;
		}
	}
	//  printf("local = %d, global = %d\n", localCount, globalCount);
	
	solver->NglobalGatherElements = globalCount;
	solver->NlocalGatherElements = localCount;
	
	if(globalCount)
		solver->o_globalGatherElementList =
		  mesh->device.malloc(globalCount*sizeof(int), solver->globalGatherElementList);
		  
	if(localCount)
		solver->o_localGatherElementList =
		  mesh->device.malloc(localCount*sizeof(int), solver->localGatherElementList);
		  
	free(localHaloFlags);
	
	return solver;
}
