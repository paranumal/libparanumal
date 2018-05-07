
// Boundary conditions
/* wall 1, outflow 2 */
#define acousticsDirichletConditions2D(bc, t, x, y, nx, ny, rM, uM, vM, rB, uB, vB) \
{                                   \
  if(bc==2){                        \
    *(rB) = -rM;                     \
    *(uB) = uM;                    \
    *(vB) = vM;                    \
  } else if(bc==1){          \
    *(rB) = rM;                     \
    *(uB) = uM - p_two*(nx*uM+ny*vM)*nx;  \
    *(vB) = vM - p_two*(nx*uM+ny*vM)*ny;  \
  }                                 \
}

