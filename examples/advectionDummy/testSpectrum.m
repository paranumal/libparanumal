
A = transpose(load('transA.dat'));
S = transpose(load('transS.dat'));
F = transpose(load('transF.dat'));
M = diag(load('M.dat'));

%% skew forms
test2 = 1
MK = 0.5*(( (M*A) - transpose(M*A) ));
MFK = 0.5*(( (M*F*A) - transpose(M*F*A) ));
MS = M*S;
test2 = 2
[dS] = eig(MS, M);
test3 = 3
[dK] = eig(MK, M);
test4 = 4
[dA] = eig(MK+MS, M);
test5 = 5
[dFA] = eig(MFK+MS, M);
test6 = 6

[min(real(dS)), max(real(dS)), min(imag(dS)), max(imag(dS))]
[min(real(dK)), max(real(dK)), min(imag(dK)), max(imag(dK))]
[min(real(dA)), max(real(dA)), min(imag(dA)), max(imag(dA))]
[min(real(dFA)), max(real(dFA)), min(imag(dFA)), max(imag(dFA))]
