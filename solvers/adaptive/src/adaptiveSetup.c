#include "adaptive.h"

/** create adaptive structure
 *
 * \param [in]  options    the options
 * \param [in]  comm       MPI communicator
 *
 * \return Initialized adaptive structure
 */
adaptive_t *adaptive_new(setupAide &options, MPI_Comm comm)
{

  printf("!!!!!!!!!!!!!   sizeof(dfloat_t) = %zu\n", sizeof(dfloat_t));
  
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //  adaptive_t *adaptive = (adaptive_t *)asd_malloc(sizeof(adaptive_t));
  adaptive_t *adaptive = new adaptive_t[1];

  adaptive->options = options;
  adaptive->allNeumann = 0; // hack for the moment (assume lambda>0)
  
  adaptive->comm = comm;
  adaptive->rank = rank;
  adaptive->size = size;
  
  occaDeviceConfig(adaptive->device, options, comm);

  adaptive->brick_n[0] = 10;
  adaptive->brick_n[1] = 10;
  adaptive->brick_n[2] = 10;
  adaptive->brick_p[0] = 1;
  adaptive->brick_p[1] = 1;
  adaptive->brick_p[2] = 1;

  options.getArgs("BOX NX", adaptive->brick_n[0]);
  options.getArgs("BOX NY", adaptive->brick_n[1]);
  options.getArgs("BOX NZ", adaptive->brick_n[2]);

  options.getArgs("BOX PERIODIC X", adaptive->brick_p[0]);
  options.getArgs("BOX PERIODIC Y", adaptive->brick_p[1]);
  options.getArgs("BOX PERIODIC Z", adaptive->brick_p[2]);

  adaptive->conn = get_connectivity(NULL, adaptive->brick_n, adaptive->brick_p);

  // This is needed for periodicity in the brick
  {
    adaptive->brick_TToC =
      (p4est_topidx_t *)asd_malloc_aligned(sizeof(p4est_topidx_t) *
					   adaptive->conn->num_trees * 3);
    const double *vertices = adaptive->conn->vertices;
    const p4est_topidx_t *tree_to_vertex = adaptive->conn->tree_to_vertex;

    for (p4est_topidx_t t = 0; t < adaptive->conn->num_trees; ++t)
      for (int d = 0; d < 3; ++d)
        adaptive->brick_TToC[3 * t + d] = (p4est_topidx_t)
          vertices[3 * tree_to_vertex[t * P4EST_CHILDREN] + d];
  }

#if HAVE_COORDINATE_MADAPTIVEING
  {
    /* let the user modify the connectivity vertices */
    ASD_ASSERT(sizeof(int) >= sizeof(p4est_topidx_t));
    for (p4est_topidx_t i = 0; i < adaptive->conn->num_vertices; i++)
      {
	double x = adaptive->conn->vertices[i * 3 + 0];
	double y = adaptive->conn->vertices[i * 3 + 1];
	double z = adaptive->conn->vertices[i * 3 + 2];
	int result =
          asd_lua_global_function_call(adaptive->prefs->L, "adaptive.conn.vertex_madaptiveing",
                                       "iddd>ddd", (int)i, x, y, z, &x, &y, &z);
	if (result != 0)
	  ASD_ABORT("Failed to move node %jd (%d)", (intmax_t)i);
	adaptive->conn->vertices[i * 3 + 0] = x;
	adaptive->conn->vertices[i * 3 + 1] = y;
	adaptive->conn->vertices[i * 3 + 2] = z;
      }
  }
#endif

  int start_level = 0;
  options.getArgs("STARTING REFINEMENT LEVEL", start_level);

  adaptive->pxest = p4est_new_ext(comm, adaptive->conn, 0, start_level, 1,
			     sizeof(quad_data_t), NULL, NULL);

  p4est_balance_ext(adaptive->pxest, P4EST_CONNECT_FULL, NULL, NULL);
  p4est_partition(adaptive->pxest, 1, NULL);
  adaptive->ghost = p4est_ghost_new(adaptive->pxest, P4EST_CONNECT_FULL);

  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  // TODO build more than one level
  adaptive->lvl = level_new(options, adaptive->pxest, adaptive->ghost, adaptive->device,
			    adaptive->brick_n, adaptive->brick_p, adaptive->brick_TToC, N, 0.01, adaptive->comm);

  occa::properties kernelInfo;
  
  kernelInfo["defines"].asObject();
  kernelInfo["defines/p_blockSize"]= KERNEL_REDUCE_LDIM;

  if(sizeof(iint_t)==4)
    kernelInfo["defines/dlong"]= "int";
  else
    kernelInfo["defines/dlong"]= "long long int";

  if(sizeof(dfloat_t)==4)
    kernelInfo["defines/dfloat"]= "float";
  else
    kernelInfo["defines/dfloat"]= "double";
  
  char fileName[BUFSIZ], kernelName[BUFSIZ];
  
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {      

      //mesh kernels
      adaptive->addScalarKernel =
        adaptive->device.buildKernel(DHOLMES "/okl/addScalar.okl",
				"addScalar",
				kernelInfo);
      
      adaptive->sumKernel =
        adaptive->device.buildKernel(DHOLMES "/okl/sum.okl",
				"sum",
				kernelInfo);

      adaptive->weightedInnerProduct1Kernel =
        adaptive->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
				"weightedInnerProduct1",
				kernelInfo);

      adaptive->weightedInnerProduct2Kernel =
        adaptive->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
				"weightedInnerProduct2",
				kernelInfo);

      adaptive->innerProductKernel =
        adaptive->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
				"innerProduct",
				kernelInfo);

      adaptive->weightedNorm2Kernel =
        adaptive->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
				"weightedNorm2",
				kernelInfo);

      adaptive->norm2Kernel =
        adaptive->device.buildKernel(DHOLMES "/okl/norm2.okl",
				"norm2",
				kernelInfo);

      adaptive->scaledAddKernel =
	adaptive->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
				"scaledAdd",
				kernelInfo);

      adaptive->dotMultiplyKernel =
	adaptive->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
				"dotMultiply",
				kernelInfo);

      adaptive->dotMultiplyAddKernel =
	adaptive->device.buildKernel(DHOLMES "/okl/dotMultiplyAdd.okl",
				"dotMultiplyAdd",
				kernelInfo);

      adaptive->dotDivideKernel =
	adaptive->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
				"dotDivide",
				kernelInfo);
    }
    MPI_Barrier(comm);
  }

  return adaptive;
}


void adaptive_free(adaptive_t *adaptive)
{
  asd_free_aligned(adaptive->brick_TToC);

  level_free(adaptive->lvl);

  p4est_ghost_destroy(adaptive->ghost);
  p4est_destroy(adaptive->pxest);
  p4est_connectivity_destroy(adaptive->conn);

  delete [] adaptive;
}
