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
  //  app_t *app = (app_t *)asd_malloc(sizeof(app_t));
  app_t *app = new app_t[1];

  occaDeviceConfig(app->device, options, comm);

  int n[DIM], p[DIM];

  n[0] = 10; n[1] = 10; n[2] = 10;
  p[0] = 1; p[1] = 1; p[2] = 1;

  options.getArgs("BOX NX", n[0]);
  options.getArgs("BOX NY", n[1]);
  options.getArgs("BOX NZ", n[2]);

  options.getArgs("BOX PERIODIC X", n[0]);
  options.getArgs("BOX PERIODIC Y", n[1]);
  options.getArgs("BOX PERIODIC Z", n[2]);

  app->conn = get_connectivity(NULL, n, p);

  return app;
}

void app_free(app_t *app)
{
  // TODO delete app->device
  // p4est_ghost_destroy(app->ghost);
  // p4est_destroy(app->pxest);
  p4est_connectivity_destroy(app->conn);

  delete [] app;

}
