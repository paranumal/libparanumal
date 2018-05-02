

// Boundary conditions
/* wall 1, inflow 2 */
#define acousticsDirichletConditions3D(bc, t, x, y, z, nx, ny, nz, rM, uM, vM, wM, rB, uB, vB, wB) \
  {									\
    if(bc==2){								\
      *(rB) = rM;							\
      *(uB) = 0.f;							\
      *(vB) = 0.f;							\
      *(wB) = 0.f;							\
    } if(bc>0){								\
      *(rB) = rM;							\
      *(uB) = uM - p_two*(nx*uM+ny*vM+nz*wM)*nx;			\
      *(vB) = vM - p_two*(nx*uM+ny*vM+nz*wM)*ny;			\
      *(wB) = wM - p_two*(nx*uM+ny*vM+nz*wM)*nz;			\
    }									\
  }
