function ellipticHexScanPlot(M, col, xlab, ylab, legloc)

forms = [{'r*-'},{'go-'},{'ks-'},{'bd-'},{'ro-'},{'c+-'},{'m-'}, {'ko-'}];
% M = load('ellipticHexScanMultigrid.results');

%M = load('ellipticHexScanJacobi.results');
maxN = max(M(:,1));

clf
hold on

for N=1:maxN
ids = find(M(:,1)==N)
[foo,ind] = sort(M(ids,2));
Msort = M(ids(ind),:);
ha = loglog(Msort(:,2), Msort(:,col), forms{N});
set(ha, 'markerfacecolor', forms{N}(1))
end

legend('N=1', 'N=2', 'N=3', 'N=4', 'N=5', 'N=6', 'N=7', 'N=8', 'N=9', 'N=10', 'location', legloc)

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
grid on
xlabel('dofs = #elements x #nodes')
ylabel(ylab)
