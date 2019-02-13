#ifndef LEVEL_H
#define LEVEL_H 1

#include "adaptive.h"

typedef struct level
{
  int N;      // order of the elements
  int Nq;     // N + 1
  int Nqk;    // dofs in the k direction (=1 for 2D, =Nq for 3D)
  int Np;     // dofs per element
  int Nfp;    // dofs per face
  int Nfaces; // faces per element

  iint_t Nmortar;     // number of mortar faces
  iint_t Ncontinuous; // number of continuous nodes
  iint_t Ncindices;   // number of continuous indices

  iint_t Kmax;        // number of elements used for the allocations
  iint_t Ktotal;      // number of total elements (these may not all be valid)
  iint_t Klocal;      // number of local elements
  iint_t Kintra;      // number of local interior (aka not mirror) elements
  iint_t Kmirror;     // number of local mirror elements
  iint_t Kuniqmirror; // number of unique local mirror elements
  iint_t Kghost;      // number of remote elements

  // p4est conn information
  occa::memory o_tree_to_vertex;
  occa::memory o_tree_vertices;

  // p4est information
  int Nn;    // number of mirror and ghost neighbors
  int *NToR; // neighbor to rank

  int8_t *EToA; // refine, none, and coarsen flag for each element

  occa::memory o_IToE;  // intra elements
  occa::memory o_MToE;  // mirror elements
  occa::memory o_UMToE; // unique mirror elements
  occa::memory o_GToE;  // ghost elements

  occa::memory o_EToL;   // element to p4est level
  occa::memory o_EToT;   // element to p4est treeid
  occa::memory o_EToX;   // element to p4est x-coord
  occa::memory o_EToY;   // element to p4est y-coord
  occa::memory o_EToZ;   // element to p4est z-coord
  occa::memory o_EToB;   // element to boundary condition
  occa::memory o_EToE;   // element to neighbor element
  occa::memory o_EToF;   // element to neighbor face
  occa::memory o_EToO;   // element to neighbor orientation
  occa::memory o_EToC;   // element to p4est lnodes face_code
  occa::memory o_EToP;   // element to periodicity mask
  occa::memory o_EToOff; // element to offset used in adaptivity

  occa::memory o_CToD_starts;  // start of indices for each continuous node
  occa::memory o_CToD_indices; // indices for continuous to discontinuous node map

  occa::memory o_MFToEM; // mortar face to minus element
  occa::memory o_MFToFM; // mortar face to minus element face
  occa::memory o_MFToEP; // mortar face to plus elements
  occa::memory o_MFToFP; // mortar face to plus elements face
  occa::memory o_MFToOP; // mortar face to plus elements orientation

  // element operators
  occa::memory o_r;
  occa::memory o_w;
  occa::memory o_D;

  // element interpolation operators
  occa::memory o_Ib;
  occa::memory o_It;

  // element L^2 projection operators
  occa::memory o_Pb;
  occa::memory o_Pt;

  // fields
  occa::memory o_q;
  occa::memory o_rhsq;

  // arrays for sending & receiving q
  occa::memory o_q_buf;

  occa::memory pin_q_recv;
  occa::memory pin_q_send;

  dfloat_t *q_recv;
  dfloat_t *q_send;

  dfloat_t *q_recv_buf;
  dfloat_t *q_send_buf;

  MPI_Request *q_recv_requests;
  MPI_Request *q_send_requests;

  MPI_Status *q_recv_statuses;
  MPI_Status *q_send_statuses;

  // geometry information
  occa::memory o_vgeo;
  occa::memory o_sgeo;
  occa::memory o_ggeo;

  // reduction buffers
  occa::memory o_red_buf[2];

  // kernels
  occa::kernel compute_Ax;// not populated yet
  occa::kernel compute_partial_Ax;
  
  occa::kernel compute_X;
  occa::kernel interp_X;
  occa::kernel coarse_X;

  occa::kernel compute_geo;
  occa::kernel compute_ics;
  occa::kernel compute_dt;
  occa::kernel compute_energy;
  occa::kernel compute_error;

  occa::kernel coarsen_fields;
  occa::kernel refine_and_fill_fields;

  occa::kernel volume_advection;
  occa::kernel mortar_advection;
  occa::kernel update_advection;
  occa::kernel zero_fields;

  occa::kernel get_mirror_fields;
  occa::kernel set_ghost_fields;

  occa::kernel reduce_min;
  occa::kernel reduce_sum;
} level_t;

level_t *level_new(setupAide &options, p4est_t *pxest,
                   p4est_ghost_t *ghost, occa::device &device,
                   int *brick_n, int *brick_p, int *brick_TToC,
                   int N, double occa_kmax_mem_frac);

void level_free(level_t *lvl);

#endif
