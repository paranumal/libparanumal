#include "adaptive.h"

/** create app structure
 *
 * \param [in]  options    the options
 * \param [in]  comm       MPI communicator
 *
 * \return Initialized app structure
 */
app_t *app_new(setupAide &options, MPI_Comm comm)
{

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  //  app_t *app = (app_t *)asd_malloc(sizeof(app_t));
  app_t *app = new app_t[1];

  occaDeviceConfig(app->device, options, comm);

  app->brick_n[0] = 10;
  app->brick_n[1] = 10;
  app->brick_n[2] = 10;
  app->brick_p[0] = 1;
  app->brick_p[1] = 1;
  app->brick_p[2] = 1;

  options.getArgs("BOX NX", app->brick_n[0]);
  options.getArgs("BOX NY", app->brick_n[1]);
  options.getArgs("BOX NZ", app->brick_n[2]);

  options.getArgs("BOX PERIODIC X", app->brick_p[0]);
  options.getArgs("BOX PERIODIC Y", app->brick_p[1]);
  options.getArgs("BOX PERIODIC Z", app->brick_p[2]);

  app->conn = get_connectivity(NULL, app->brick_n, app->brick_p);

  // This is needed for periodicity in the brick
  {
    app->brick_TToC =
      (p4est_topidx_t *)asd_malloc_aligned(sizeof(p4est_topidx_t) *
					   app->conn->num_trees * 3);
    const double *vertices = app->conn->vertices;
    const p4est_topidx_t *tree_to_vertex = app->conn->tree_to_vertex;

    for (p4est_topidx_t t = 0; t < app->conn->num_trees; ++t)
      for (int d = 0; d < 3; ++d)
        app->brick_TToC[3 * t + d] = (p4est_topidx_t)
          vertices[3 * tree_to_vertex[t * P4EST_CHILDREN] + d];
  }

#if HAVE_COORDINATE_MAPPING
  {
    /* let the user modify the connectivity vertices */
    ASD_ASSERT(sizeof(int) >= sizeof(p4est_topidx_t));
    for (p4est_topidx_t i = 0; i < app->conn->num_vertices; i++)
      {
	double x = app->conn->vertices[i * 3 + 0];
	double y = app->conn->vertices[i * 3 + 1];
	double z = app->conn->vertices[i * 3 + 2];
	int result =
          asd_lua_global_function_call(app->prefs->L, "app.conn.vertex_mapping",
                                       "iddd>ddd", (int)i, x, y, z, &x, &y, &z);
	if (result != 0)
	  ASD_ABORT("Failed to move node %jd (%d)", (intmax_t)i);
	app->conn->vertices[i * 3 + 0] = x;
	app->conn->vertices[i * 3 + 1] = y;
	app->conn->vertices[i * 3 + 2] = z;
      }
  }
#endif

  int start_level = 0;
  options.getArgs("STARTING REFINEMENT LEVEL", start_level);

  app->pxest = p4est_new_ext(comm, app->conn, 0, start_level, 1,
			     sizeof(quad_data_t), NULL, NULL);

  p4est_balance_ext(app->pxest, P4EST_CONNECT_FULL, NULL, NULL);
  p4est_partition(app->pxest, 1, NULL);
  app->ghost = p4est_ghost_new(app->pxest, P4EST_CONNECT_FULL);

  int N;
  options.getArgs("POLYNOMIAL DEGREE", N);
  // TODO build more than one level
  app->lvl = level_new(options, app->pxest, app->ghost, app->device,
      app->brick_n, app->brick_p, app->brick_TToC, N);

  int blockSize = 256;
  occa::properties kernelInfo;
  
  kernelInfo["defines"].asObject();
  kernelInfo["defines/p_blockSize"]= blockSize;

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
      app->addScalarKernel =
        app->device.buildKernel(DHOLMES "/okl/addScalar.okl",
				"addScalar",
				kernelInfo);
      
      app->sumKernel =
        app->device.buildKernel(DHOLMES "/okl/sum.okl",
				"sum",
				kernelInfo);

      app->weightedInnerProduct1Kernel =
        app->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
				"weightedInnerProduct1",
				kernelInfo);

      app->weightedInnerProduct2Kernel =
        app->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
				"weightedInnerProduct2",
				kernelInfo);

      app->innerProductKernel =
        app->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
				"innerProduct",
				kernelInfo);

      app->weightedNorm2Kernel =
        app->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
				"weightedNorm2",
				kernelInfo);

      app->norm2Kernel =
        app->device.buildKernel(DHOLMES "/okl/norm2.okl",
				"norm2",
				kernelInfo);

      app->scaledAddKernel =
	app->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
				"scaledAdd",
				kernelInfo);

      app->dotMultiplyKernel =
	app->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
				"dotMultiply",
				kernelInfo);

      app->dotMultiplyAddKernel =
	app->device.buildKernel(DHOLMES "/okl/dotMultiplyAdd.okl",
				"dotMultiplyAdd",
				kernelInfo);

      app->dotDivideKernel =
	app->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
				"dotDivide",
				kernelInfo);
    }
    MPI_Barrier(comm);
  }  
  return app;
}


void app_free(app_t *app)
{
  asd_free_aligned(app->brick_TToC);

  level_free(app->lvl);

  p4est_ghost_destroy(app->ghost);
  p4est_destroy(app->pxest);
  p4est_connectivity_destroy(app->conn);

  delete [] app;
}
