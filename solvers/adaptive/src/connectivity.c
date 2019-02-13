#include "adaptive.h"

// {{{ Connectivity
/** Get a new connectivity
 *
 * \param [in]  name    p4est name or file name of connectivity;
 *                      or NULL for brick connectivity
 * \param [in]  n[DIM]  number (per dimension) of quads
 * \param [in]  p[DIM]  Boolean (per dimension) indicating periodicity of brick
 *
 * \return Initialized app structure
 *
 */
p4est_connectivity_t *get_connectivity(const char *name, int *n, int *p)
{
  p4est_connectivity_t *conn = NULL;

  if (name)
  {
    FILE *stream = fopen(name, "r");
    if (stream)
    {
      if (fclose(stream))
        ASD_ABORT("Failed fclose on: %s", name);

      conn = p4est_connectivity_read_inp(name);
    }
    else
    {
      conn = p4est_connectivity_new_byname(name);
    }

    if (conn == NULL)
      ASD_ABORT("Failed get valid pxest connectivity \"%s\"", name);
  }
  else
  {
    conn = p4est_connectivity_new_brick(n[0], n[1],
#if DIM == 3
                                        n[2],
#endif
                                        p[0], p[1]
#if DIM == 3
                                        ,
                                        p[2]
#endif
                                        );
  }

  return conn;
}
// }}}
