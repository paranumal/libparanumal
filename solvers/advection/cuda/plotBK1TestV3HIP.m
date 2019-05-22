
data = load('results/testMassMatrixMultiplyVT_RADEONVII_HIP_V3.dat')

dNq = 1;

ids = find(data(:,1)+dNq==data(:,2));

data = data(ids,:);

clf
hold on
maxBW = 0;
for mode=1:6
  ids = find(data(:,end)==mode)
  scatter(data(ids,1)-1, data(ids,8), 'filled')
  maxBW = max(maxBW, data(ids,9));
end

plot(data(ids,1)-1, maxBW, 'r-', 'LineWidth', 2)

hold off

axis([1, 11, 0 800])

  
grid on
xlabel('Element Degree', 'FontSize', 14)
ylabel('Estimated Bandwidth (GB/s)', 'FontSize', 14)
title('BK1:HIP:RADEON VII:~8M DOFS')
ha = legend('Odd-even + OP in registers', ...
	    'Odd-even + OP in constant cache', ...
	    'Odd-even + OP in shared cache', ...
	    'Odd-even + OP in global', ...
	    'Monolithic + OP in Global', ...
	    'Monolithic + OP in constant cache', ...
	    'hipMemcpy BW', 'location', 'southwest');
set(ha, 'FontSize', 12)

print('-dpdf', 'testMassMatrixMultiplyVT_RADEONVII_HIP_V3.pdf', '-bestfit')
