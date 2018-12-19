

 TitanV = load('bandwidthTesterTitanV.dat');

 figure(1)
 plot(TitanV(:,1), TitanV(:,3))
 hold on
 maxBW = max(TitanV(:,3));
 plot([TitanV(1,1),TitanV(end,1)], maxBW*[.8,.8], 'r-')
 hold off
 xlabel('Bytes')
 ylabel('GB/s');
 grid minor
 print('-dpdf', 'bandwidthTesterBytesTitanV.pdf');

 figure(2)
 plot(TitanV(:,1)/64, TitanV(:,3))
 hold on
 plot([TitanV(1,1),TitanV(end,1)]/64, maxBW*[.8,.8], 'r-')
 hold off
 xlabel('SEM Nodes (64 bytes for BK5)')
 ylabel('GB/s');
 grid minor
 print('-dpdf', 'bandwidthTesterNodesTitanV.pdf');

axis([0 2e5 0 600]) 
 print('-dpdf', 'bandwidthTesterNodesZoomTitanV.pdf');

figure(3)

semilogy(TitanV(:,1)/64, TitanV(:,2)./(TitanV(:,1)/64), 'linewidth', 1.8)

minTimePerNode = min(TitanV(:,2)./(TitanV(:,1)/64));

 hold on
 plot([TitanV(1,1),TitanV(end,1)]/64, minTimePerNode*[1./.8,1./.8], 'r-')
 hold off
 xlabel('SEM Nodes (64 bytes for BK5)')
 ylabel('copy time/node (s)');
 grid minor
 set(gca, 'xlim', [0 2000000])
   set(gca, 'ylim', [1e-10 3e-10])
 print('-dpdf', 'bandwidthTesterNodesCopyTimeTitanV.pdf');

%axis([0 2e5 0 600]) 
% print('-dpdf', 'bandwidthTesterNodesZoomTitanV.pdf');
