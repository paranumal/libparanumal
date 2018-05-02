
// Boundary conditions
/* wall 1, outflow 2 */
#define acousticsDirichletConditions2D(bc, t, x, y, nx, ny, rM, uM, vM, rB, uB, vB) \
{                                   \
  if(bc==2){                        \
    *(rB) = rM;                     \
    *(uB) = 0.f;                    \
    *(vB) = 0.f;                    \
  } else if(bc==1){          \
    *(rB) = rM;                     \
    *(uB) = uM - (nx*uM+ny*vM)*nx;  \
    *(vB) = vM - (nx*uM+ny*vM)*ny;  \
  }                                 \
}

