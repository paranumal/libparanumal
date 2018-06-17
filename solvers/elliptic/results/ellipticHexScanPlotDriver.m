
J = load('ellipticHexScanJacobi.results');

M = load('ellipticHexScanMultigrid.results');

figure(1)
ellipticHexScanPlot(J, 5, 'dofs', '(solve time)/dof', 'northeast')
title('Hexes: Jacobi Preconditioner');
myprint('ellipticPoissonSolverHexJacobiSolveTime.pdf')

figure(2)
ellipticHexScanPlot(M, 5, 'dofs', '(solve time)/dof', 'northeast')
title('Hexes: Hybrid AMG/pMG Preconditioner');
myprint('ellipticPoissonSolverHexHybridMultigridSolveTime.pdf')

figure(3)
ellipticHexScanPlot(J, 6, 'dofs', '(iterations x dofs)/(solve time)', 'southeast')
title('Hexes: Jacobi Preconditioner');
myprint('ellipticPoissonSolverHexHybridJacobiThroughput.pdf')

figure(4)
ellipticHexScanPlot(M, 6, 'dofs', '(iterations x dofs)/(solve time)', 'southeast')
title('Hexes: Hybrid AMG/pMG Preconditioner');
myprint('ellipticPoissonSolverHexHybridMultigridThroughput.pdf')

if(0)
J = load('ellipticTetScanJacobi.results');

M = load('ellipticTetScanMultigrid.results');


figure(6)
ellipticHexScanPlot(M, 5, 'dofs', '(solve time)/dof', 'southeast')
title('Tets: Hybrid AMG/pMG Preconditioner');
myprint('ellipticPoissonSolverTetMultigridSolveTime.pdf')

figure(8)
ellipticHexScanPlot(M, 6, 'dofs', '(iterations x dofs)/(solve time)', 'southeast')
title('Tets: Hybrid AMG/pMG Preconditioner');
myprint('ellipticPoissonSolverTetHybridMultigridThroughput.pdf')

figure(5)
ellipticHexScanPlot(J, 5, 'dofs', '(solve time)/dof', 'southeast')
title('Tets: Jacobi Preconditioner');
myprint('ellipticPoissonSolverTetJacobiSolveTime.pdf')

figure(7)
ellipticHexScanPlot(J, 6, 'dofs', '(iterations x dofs)/(solve time)', 'southeast')
title('Tets: Jacobi Preconditioner');
myprint('ellipticPoissonSolverTetJacobiThroughput.pdf')

end
