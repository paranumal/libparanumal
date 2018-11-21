
function rooflinePlot(kernelName)

  perf = load(kernelName);

  clf

  scatter(perf(:,1),perf(:,2), 'filled', 'ko');

  hold on;

  set(gca, 'FontSize', 14);
  xlabel('Arithmetic Intensity (FP64 flops)/deviceMemoryBytes', 'FontSize', 14);
  ylabel('FP64 GFLOPS/s', 'FontSize', 14);
  title(kernelName, 'FontSize', 16, 'Interpreter', 'None');
  
  % superimpose roofline
  hold on;
  plot(perf(:,1),perf(:,1).*perf(:,3),'r*', 'LineWidth', 2);
  axis([0,ceil(max(perf(:,1))),0,max(max(perf(:,1).*perf(:,3),perf(:,2)))]);
