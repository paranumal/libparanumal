#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

/** Get a new p4est connectivity
 *
 * \param [in]  name    p4est name or file name of connectivity;
 *                      or NULL for brick connectivity
 * \param [in]  n[DIM]  number (per dimension) of quads
 * \param [in]  p[DIM]  Boolean (per dimension) indicating periodicity of brick
 *
 * \return Initialized app structure
 *
 */
p4est_connectivity_t *get_connectivity(const char *name, int *n, int *p);

#endif
