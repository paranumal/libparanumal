cudaFlags = -arch=sm_70   --use_fast_math -Xptxas -v,-dlcm=cg,-O3 -Xcompiler -O3
#-g -G

advectionInvertMassMatrix: advectionInvertMassMatrix.cu
	nvcc -o advectionInvertMassMatrix $(cudaFlags) advectionInvertMassMatrix.cu

clean:
	rm -f advectionInvertMassMatrix
