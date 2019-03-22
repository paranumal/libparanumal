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


#include <stdio.h>
#include <stdlib.h>

#include "occa.hpp"

#define int int
#define intString "int"

#if 1
#define dfloat double
#define dfloatString "double"
#else
#define dfloat float
#define dfloatString "float"
#endif

// this is for optimizing ellpack matrix format

#define maxNNZperROW 25


int main(int argc, char **argv){

	int Ntests = 10;
	// N regulates the size of the matrix. The matrix is random.
	
	int  N = atoi(argv[1]);
	
	int seed = time(NULL);
	srand(seed);
	//size of register array (in doubles) per thread
	//int reg = atoi(argv[2]);
	
	
	// this matrix format has two arrays
	// maxNNZperROW*N values (double) and maxNNZperROW*N column indices (int)
	double *h_val = (double*) calloc(maxNNZperROW*N , sizeof(double));
	int    *h_ind = (int*) calloc(maxNNZperROW*N , sizeof(int));
	double *h_vector = (double*) calloc(N, sizeof(double));
	double *h_CPUresult = (double*) calloc(N, sizeof(double));
	for(int n=0;n<N*maxNNZperROW;++n){
		h_ind[n] =0;
	}
	
	// put some random stuff in the arrays
	int Nd = N/5000;
	int currentStart = 0;
	int m;
	for(int n=0;n<N;++n){
	
		if ((n%Nd==0)&&(n!=0)) {
			//	printf("adding %d to %d, n= %d \n",Nd, currentStart, n);
			currentStart +=Nd;
			
		}
		h_vector[n] = (double) rand()/ (RAND_MAX-8);
		int it = 0;
		int lastindex = -1;
		int howmany = rand()%maxNNZperROW;
		if (howmany ==0) howmany++;
		//	printf("this is row %d, maxx nnz %d \n", n, howmany);
		m=0;
		while (m<howmany){
		
			double val =  (double)rand() / RAND_MAX;
			int index = rand()%Nd;
			//		printf("got %d \n", index);
			//		printf("Nd  %d row %d currentStart %d \n", Nd, n, currentStart);
			//index must be bigger than the last index but smaller than maxNNZ-m;
			while ((index<=lastindex))
			{
				//	printf()
				index = rand()%Nd;
				//		printf("got (repeat) %d \n", index);
			}
			if (index >=Nd-(howmany)+m)
			{
				index = Nd-(howmany)+m;
				
				for (int i = index; i<Nd; i++)
				{
					double val =  (double)rand() / RAND_MAX;
					h_val[n*maxNNZperROW+m] = val;
					h_ind[n*maxNNZperROW+m] = currentStart+i+1;
					
					m++;
					
				}
				m=N;
				
			}
			else{
			
				h_val[n*maxNNZperROW+m] = val;
				h_ind[n*maxNNZperROW+m] = currentStart+index+1;
				
				lastindex=index;
				
				
				m++;
				
				
			}
			
		}
		
		
	}
	
	
	/*	for(int n=0;n<N;++n){
		printf("row %d: ", n);
		for(int m=0;m<maxNNZperROW;++m){
			printf(" %d ", h_ind[n*maxNNZperROW+m]);
		}
		printf("\n");
}

				for(int n=0;n<N;++n){
					printf("%f ", h_vector[n]);
			
}*/
	double normsum   =0.0;
	for (int n=0; n<N; ++n)
	{
		//row by row
		double rowsum = 0.0;
		for (int m=0; m<maxNNZperROW; ++m)
		{
			int colindex = h_ind[n*maxNNZperROW+m];
			if (colindex !=0)
			{
				rowsum += h_val[m+n*maxNNZperROW]*h_vector[colindex-1];
				//	printf("multiplying nnz in row %d column %d by vec[%d] \n", n, colindex,colindex  );
			}
		}
		//	h_CPUresult[n]=rowsum;
		normsum+=(rowsum*rowsum);
	}
	printf("vector norm %f\n", sqrt(normsum));
	double *h_out = (double*) calloc(N, sizeof(double));
	printf("chujemuje \n");
	
	
	
	char deviceConfig[BUFSIZ];
	sprintf(deviceConfig, "mode = CUDA, deviceID = 0");
	
	occa::device device;
	device.setup(deviceConfig);
	occa::memory o_val = device.malloc(N*maxNNZperROW*sizeof(double), h_val);
	occa::memory o_vector = device.malloc(N*sizeof(double), h_vector);
	occa::memory o_ind = device.malloc(maxNNZperROW*N*sizeof(int), h_ind);
	occa::memory o_out = device.malloc(N*maxNNZperROW*sizeof(double), h_out);
	printf("chujemuje 2\n");
	int blockSize = 128;
	occa::kernelInfo kernelInfo;
	kernelInfo.addDefine("p_rowSize", maxNNZperROW);
	//kernelInfo.addDefine("p_blockSize", blockSize);
	//kernelInfo.addDefine("p_Nblocks", (int) N/blockSize+1);
	//for v6
	printf("number of blocks %d blockSize %d \n", (int)N/Nd, Nd);
	kernelInfo.addDefine("p_blockSize", Nd);
	kernelInfo.addDefine("p_Nblocks", N/Nd);
	
	//	kernelInfo.addDefine("p_Nq", N+1);
	//	kernelInfo.addDefine("p_Np", (N+1)*(N+1)*(N+1));
	//	kernelInfo.addDefine("p_gjNp", (Nq+1)*(Nq+1)*(Nq+1));
	
	//kernelInfo.addDefine("elPerBlock", Np*2+7*gjNp);
	
	kernelInfo.addParserFlag("automate-add-barriers", "disabled");
	kernelInfo.addCompilerFlag("  --compiler-options -O3");
	printf("chujemuje \n");
	occa::kernel testMatvec
	= device.buildKernelFromSource(LIBP_DIR "/okl/testHYBmatvec.okl",
"testHYBmatvec_v6",
	                               kernelInfo);
	                               
	/*int N,
	                            double * values,
	                            int * indices,
	                            double * vector,
	                            double * result*/
	occa::streamTag startTag = device.tagStream();
	
	for(int test=0;test<Ntests;++test){
		testMatvec(N, o_val, o_ind, o_vector, o_out);
	}
	
	occa::streamTag stopTag = device.tagStream();
	o_out.copyTo(h_out);
	double GPUnorm = 0.0f;
	for (int n=0; n<N; ++n)
	{
		GPUnorm += (h_out[n]*h_out[n]);
		//	if (h_out[n]!=h_CPUresult[n])
		//	printf("GPU=%d in %f CPU %f \n", n, h_out[n], h_CPUresult[n]);
	}
	printf("CPU vector norm %15.15f GPU vector norm %15.15f difference = %15.15E\n", sqrt(GPUnorm), sqrt(normsum), sqrt(GPUnorm)- sqrt(normsum));
	//double copyElapsed = device.timeBetween(startTag, stopTag);
	
	//Nbytes =(sizeof(double)/2)*(Np*2 +7*gjNp);
	//double copyBandwidth = Nelements*((Nbytes*Ntests*2.)/(1024.*1024.*1024.*copyElapsed));
	
	//printf("time %8.8f bw %17.15E \n", copyElapsed, copyBandwidth);
	
}
