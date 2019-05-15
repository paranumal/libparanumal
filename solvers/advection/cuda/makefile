# -cg 42%
cudaFlags = -Dcomp_Nq=$(comp_Nq) -Dcomp_cubNq=$(comp_cubNq) -Dcomp_MAX_ITERATIONS=$(comp_MAX_ITERATIONS) -arch=sm_70   --use_fast_math -Xptxas -v,-dlcm=ca,-O3 -Xcompiler -O3
#-g -G

advectionInvertMassMatrix: advectionInvertMassMatrixV2.cu
	nvcc -o advectionInvertMassMatrix $(cudaFlags) advectionInvertMassMatrixV2.cu

advectionMassMatrixMultiply: advectionMassMatrixMultiply.cu
	nvcc -o advectionMassMatrixMultiply $(cudaFlags) advectionMassMatrixMultiply.cu

massMatrixMultiplyVT: massMatrixMultiplyVT.cu
	nvcc -o massMatrixMultiplyVT $(cudaFlags) massMatrixMultiplyVT.cu


roofline: roofline.cu
	nvcc -o roofline $(cudaFlags) roofline.cu


clean:
	rm -f advectionInvertMassMatrix roofline
