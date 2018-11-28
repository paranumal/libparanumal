
function rooflinePlot(kernelName)

  ['!' kernelName '!']
  perf = load(kernelName)

  clf
  
  ha =   scatter(perf(:,1),perf(:,2), 'filled', 'markeredgecolor', 'r', 'markerfacecolor', 'k');

  hold on;

  set(gca, 'FontSize', 14);
  xlabel('Arithmetic Intensity (FP64 flops)/deviceMemoryBytes', 'FontSize', 14);
  ylabel('FP64 GFLOPS/s', 'FontSize', 14);
title(kernelName, 'FontSize', 16, 'Interpreter', 'None');
  
  % superimpose roofline
  [foo,ids] = sort(perf(:,1), 'ascend');
  perf(ids,:)

  maxBW = max(perf(:,end))
  
  estPerf = min(perf(ids,4));
%%TW  plot(perf(ids,1),min(perf(ids,1).*perf(ids,3),estPerf),'r*-', 'LineWidth', 2);

  maxAI = max(perf(:,1));
  plot([0,maxAI], [0,maxBW*maxAI], 'b-', 'LineWidth', 2)
  
%  plot(perf(ids,1),ones(size(ids))*estPerf,'r-', 'LineWidth', 2);
%%  axis([0,ceil(max(perf(:,1))),0,max(max(perf(:,1).*perf(:,3),perf(:,2)))]);
  hold off;
  grid minor
  
  pause(1);
