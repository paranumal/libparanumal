
Ax = load('ellipticHexScanBK5V100.results');

figure(1)
ellipticHexScanPlot(Ax, 5, 'dofs', '(Ax time)/dof', 'northeast')
title('Hexes: Ax time per dof');
myprint('ellipticHexBK5AxTime.pdf')

figure(2)
ellipticHexScanPlot(Ax, 6, 'dofs', '(dofs)/(Ax time)', 'southeast')
title('Hexes: Ax throughput ');
myprint('ellipticHexBK5AxThroughput.pdf')


