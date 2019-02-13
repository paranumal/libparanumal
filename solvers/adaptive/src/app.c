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
  app_t *app = (app_t *)asd_malloc(sizeof(app_t));

  occaDeviceConfig(app->device, options, comm);

  return app;
}

void app_free(app_t *app)
{
  // TODO delete app->device
}
