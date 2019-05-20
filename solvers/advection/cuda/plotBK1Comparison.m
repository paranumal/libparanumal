
data = load('BK1TitanVComparison.dat');

dNq = 1;

ids = find(data(:,1)+dNq==data(:,2));

data = data(ids,:);

clf
hold on
maxBW = 0;
for mode=1:3
  ids = find(data(:,end)==mode)
  scatter(data(ids,1), data(ids,8), 'filled')
  maxBW = max(maxBW, data(ids,9));
end

plot(data(ids,1), maxBW, 'b-')

plot([2,12], [653,653], 'r--', 'LineWidth', 2)


hold off

axis([2, 12, 0 700])

  
grid on
xlabel('#GLL', 'FontSize', 14)
ylabel('Estimated Bandwidth (GB/s)', 'FontSize', 14)
title('BK1:CUDA:Titan V:~4M DOFS')
ha = legend('Odd-even + matrix in registers', 'Odd-even + matrix in constant cache', 'Matrix in constant cache', 'cudaMemcpy BW', 'Manufacturer peak BW', 'location', 'southeast');
set(ha, 'FontSize', 14)
