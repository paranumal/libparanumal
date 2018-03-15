
function writeNodeData3D(inN)

Globals3D;
N = inN;
[r,s,t] = Nodes3D(N);
[r,s,t] = xyztorst(r,s,t);

Np = length(r);
Nfp = (N+1)*(N+2)/2;
Nfaces = 4;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(t+1) < NODETOL)';
faceNodes2   = find( abs(s+1) < NODETOL)';
faceNodes3   = find( abs(r+s+t+1) < NODETOL)';
faceNodes4   = find( abs(r+1) < NODETOL)';
faceNodes  = [faceNodes1;faceNodes2;faceNodes3;faceNodes4]';

V = Vandermonde3D(N, r, s, t);
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
LIFT = Lift3D(N, faceNodes, r, s, t);

fname = sprintf('tetN%02d.dat', N);

fid = fopen(fname, 'w');

fprintf(fid, '%% degree N\n');
fprintf(fid, '%d\n', N);
fprintf(fid, '%% number of nodes\n');
fprintf(fid, '%d\n', Np);
fprintf(fid, '%% node coordinates\n');
for n=1:Np
    fprintf(fid, '%17.15E %17.15E %17.15E\n', r(n), s(n), t(n));
end

fprintf(fid, '%% r collocation differentation matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', Dr(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% s collocation differentation matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', Ds(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% t collocation differentation matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', Dt(n,m));
    end
    fprintf(fid, '\n');
end

MM = inv(transpose(V))/V;
fprintf(fid, '%% reference mass matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', MM(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% faceNodes\n');
for f=1:Nfaces
    for m=1:Nfp
        fprintf(fid, '%d ', faceNodes(m,f)-1); %% adjust for 0-indexing
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% LIFT matrix\n');
for n=1:Np
    for m=1:Nfp*Nfaces
        fprintf(fid, '%17.15E ', LIFT(n,m));
    end
    fprintf(fid, '\n');
end

%% compute equispaced nodes on equilateral triangle
[plotR,plotS,plotT] = EquiNodes3D(N+4);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
%plotEToV = delaunay3(plotR,plotS,plotT)-1;
plotEToV = delaunayFixVolume(plotR,plotS,plotT)-1;

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1);

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = Vandermonde3D(N, plotR,plotS,plotT)/V;

%% output plot nodes
fprintf(fid, '%% number of plot nodes\n');
fprintf(fid, '%d\n', plotNp);
fprintf(fid, '%% plot node coordinates\n');
for n=1:plotNp
    fprintf(fid, '%17.15E %17.15E %17.15E\n', plotR(n), plotS(n), plotT(n));
end

%% output plot interpolation matrix
fprintf(fid, '%% plot node interpolation matrix\n');
for n=1:plotNp
    for m=1:Np
        fprintf(fid, '%17.15E ', plotInterp(n,m));
    end
    fprintf(fid, '\n');
end

%% output plot triangulation
fprintf(fid, '%% number of plot elements\n');
fprintf(fid, '%d\n', plotNelements);

fprintf(fid, '%% number of vertices per plot elements\n');
fprintf(fid, '%d\n', size(plotEToV,2));

fprintf(fid, '%% triangulation of plot nodes\n');
for n=1:plotNelements
    fprintf(fid, '%d %d %d %d\n' ,...
        plotEToV(n,1),plotEToV(n,2),plotEToV(n,3),plotEToV(n,4));
end

%% output cubature arrays

if N < 7
    [cubr cubs cubt cubw] = tet_cubature(2*N+1);
    Vq = Vandermonde3D(N,cubr,cubs,cubt)/V;
    Pq = V*V'*Vq'*diag(cubw);
    
    cV = Vandermonde3D(N, cubr, cubs, cubt);
    [cVr,cVs,cVt] = GradVandermonde3D(N, cubr, cubs, cubt);
    cubDrT = V*transpose(cVr)*diag(cubw);
    cubDsT = V*transpose(cVs)*diag(cubw);
    cubDtT = V*transpose(cVt)*diag(cubw);
    
    Ncub = length(cubw);
    fprintf(fid, '%% number of volume cubature nodes\n');
    fprintf(fid, '%d\n', Ncub);
    
    fprintf(fid, '%% cubature node coordinates\n');
    for n=1:Ncub
        fprintf(fid, '%17.15E %17.15E %17.15E %17.15E\n', cubr(n), cubs(n), cubt(n), cubw(n));
    end
    
    fprintf(fid,'%% cubature interpolation matrix\n');
    for i=1:Ncub
        for j = 1:Np
            fprintf(fid, '%17.15E ', Vq(i,j));
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '%% cubature r weak differentiation matrix\n');
    for n=1:Np
        for m=1:Ncub
            fprintf(fid, '%17.15E ', cubDrT(n,m));
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '%% cubature s weak differentiation matrix\n');
    for n=1:Np
        for m=1:Ncub
            fprintf(fid, '%17.15E ', cubDsT(n,m));
        end
        fprintf(fid, '\n');
    end
    
    
    fprintf(fid, '%% cubature t weak differentiation matrix\n');
    for n=1:Np
        for m=1:Ncub
            fprintf(fid, '%17.15E ', cubDtT(n,m));
        end
        fprintf(fid, '\n');
    end
    
    
    fprintf(fid,'%% cubature projection matrix\n');
    for i=1:Np
        for j = 1:Ncub
            fprintf(fid, '%17.15E ', Pq(i,j));
        end
        fprintf(fid, '\n');
    end
    
    
    
    % Surface Cubature
    [cubx,cuby,cubw] = Cubature2D(3*N); 
    Nfi = length(cubx);
    ir = [cubx,         cubx,        cubx,                  -ones(Nfi,1)];
    is = [cuby,        -ones(Nfi,1), cuby,                   cubx];
    it = [-ones(Nfi,1), cuby,       -(ones(Nfi,1)+cubx+cuby),cuby];
    iw = [cubw,cubw,cubw,cubw];
    
    sV = Vandermonde3D(N, ir(:), is(:),it(:));
    sInterp = sV/V;
    
    interp = [sInterp(1:Nfi,faceNodes(:,1));...
        sInterp(Nfi+1:2*Nfi,faceNodes(:,2));...
        sInterp(2*Nfi+1:3*Nfi,faceNodes(:,3));...
        sInterp(3*Nfi+1:4*Nfi,faceNodes(:,4))];
    fprintf(fid, '%% number of surface integration nodes per face\n');
    fprintf(fid, '%d\n', length(cubx));
    fprintf(fid, '%% surface integration interpolation matrix\n');
    for n=1:Nfi*Nfaces
        for m=1:Nfp
            fprintf(fid, '%17.15E ', interp(n,m));
        end
        fprintf(fid, '\n');
    end
    
    bInterp = [];
    bInterp(1:Nfi,1:Nfp) = interp(1:Nfi,:);
    bInterp(Nfi+1:2*Nfi,Nfp+1:2*Nfp) = interp(Nfi+1:2*Nfi,:);
    bInterp(2*Nfi+1:3*Nfi,2*Nfp+1:3*Nfp) = interp(2*Nfi+1:3*Nfi,:);
    bInterp(3*Nfi+1:4*Nfi,3*Nfp+1:4*Nfp) = interp(3*Nfi+1:4*Nfi,:);
    
    % integration node lift matrix
    iLIFT = V*V'*sInterp'*diag(iw(:));
    size(iLIFT)
    size(interp)
    max(max(abs(iLIFT*bInterp-LIFT)))
    
    fprintf(fid, '%% surface integration lift matrix\n');
    for n=1:Np
        for m=1:Nfi*Nfaces
            fprintf(fid, '%17.15E ', iLIFT(n,m));
        end
        fprintf(fid, '\n');
    end
    
end
%% BB

addpath('./bern')
tol=1e-6;

[r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
Np = (N+1)*(N+2)*(N+3)/6;
Nfaces = 4;
Nfp = (N+1)*(N+2)/2;

V = Vandermonde3D(N,r,s,t);
[VB Vr Vs Vt V1 V2 V3 V4] = bern_basis_tet(N,r,s,t);
VB(abs(VB)<tol) = 0;
invVB = inv(VB);

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
VB2D = bern_basis_tri(N,r2D,s2D);

%% write VDM for conversion


fprintf(fid, '%% Bernstein VDM matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', VB(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% inverse Bernstein VDM matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', invVB(n,m));
    end
    fprintf(fid, '\n');
end

%% write Dmats

D1 = VB\V1; D1(abs(D1)< 0.5) = 0; %all the entries are integers
D2 = VB\V2; D2(abs(D2)< 0.5) = 0;
D3 = VB\V3; D3(abs(D3)< 0.5) = 0;
D4 = VB\V4; D4(abs(D4)< 0.5) = 0;

D1ids = zeros(Np,4);
D2ids = zeros(Np,4);
D3ids = zeros(Np,4);
D4ids = zeros(Np,4);
Dvals = zeros(Np,4);
D2vals = zeros(Np,4);
D3vals = zeros(Np,4);
D4vals = zeros(Np,4);
for i = 1:Np
    tmp = find(D1(i,:));
    Dvals(i,1:length(tmp)) = D1(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D1ids(i,:) = tmp;
    
    
    tmp = find(D2(i,:));
    D2vals(i,1:length(tmp)) = D2(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D2ids(i,:) = tmp;
    
    tmp = find(D3(i,:));
    D3vals(i,1:length(tmp)) = D3(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D3ids(i,:) = tmp;
    
    tmp = find(D4(i,:));
    D4vals(i,1:length(tmp)) = D4(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D4ids(i,:) = tmp;
end

fprintf(fid, '%% sparse BB D1 matrix ids\n');
for n=1:Np
    for m=1:4
        fprintf(fid, '%d ', D1ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D2 matrix ids\n');
for n=1:Np
    for m=1:4
        fprintf(fid, '%d ', D2ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D3 matrix ids\n');
for n=1:Np
    for m=1:4
        fprintf(fid, '%d ', D3ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D4 matrix ids\n');
for n=1:Np
    for m=1:4
        fprintf(fid, '%d ', D4ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D1,D2,D3,D4 matrix values\n');
for n=1:Np
    for m=1:4
        fprintf(fid, '%17.15E ', Dvals(n,m));
    end
    fprintf(fid, '\n');
end

%% LIFT

fmask1   = find( abs(t+1) < tol)';
fmask2   = find( abs(s+1) < tol)';
fmask3   = find( abs(r+s+t+1) < tol)';
fmask4   = find( abs(r+1) < tol)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

Emat = zeros(Np, Nfaces*Nfp);
% face 1
faceR = r(Fmask(:,1));
faceS = s(Fmask(:,1));
V2D = Vandermonde2D(N, faceR, faceS);
massFace1 = inv(V2D*V2D');
Emat(Fmask(:,1),1:Nfp) = massFace1;
% face 2
faceR = r(Fmask(:,2));
faceT = t(Fmask(:,2));
V2D = Vandermonde2D(N, faceR, faceT);
massFace2 = inv(V2D*V2D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massFace2;
% face 3
faceR = r(Fmask(:,3));
faceS = s(Fmask(:,3));
V2D = Vandermonde2D(N, faceR, faceS);
massFace3 = inv(V2D*V2D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massFace3;
% face 4
faceS = s(Fmask(:,4));
faceT = t(Fmask(:,4));
V2D = Vandermonde2D(N, faceS, faceT);
massFace4 = inv(V2D*V2D');
Emat(Fmask(:,4),3*Nfp+1:4*Nfp) = massFace4;
% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
% convert to BB
LIFT = VB\(LIFT * kron(eye(Nfaces),VB2D));

EL1 = [];
for j = 0:N
    l(j+1) = (-1)^j * nchoosek(N,j)/(1+j);
    Ej = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-j,r2D,s2D);
    EL1 = [EL1; l(j+1)*Ej'];
end
EL1(abs(EL1)<tol) = 0;


[r2Dp1 s2Dp1] = Nodes2D(N+1); [r2Dp1 s2Dp1] = xytors(r2Dp1,s2Dp1);
ENp1 = bern_basis_tri(N+1,r2Dp1,s2Dp1)\bern_basis_tri(N,r2Dp1,s2Dp1);

L0 = (N+1)^2/2 * ENp1'*ENp1;
L0(abs(L0)<tol) = 0;

EL = zeros(Np,Nfaces*Nfp);
EL(:,1:Nfp) = EL1;
for f = 2:Nfaces
    ids = (1:Nfp)+(f-1)*Nfp;
    diff = sum(abs( - LIFT(:,ids)),2);
    p = zeros(Np,1);
    for i = 1:Np
        iid = find(sum(abs(repmat(LIFT(i,1:Nfp),Np,1)-LIFT(:,ids)),2)<1e-8);
        p(iid) = i;
    end
    EL(:,ids) = EL1(p,:);
end

EL(abs(EL)<tol) = 0;

% accuracy check
%if norm(LIFT-EL*kron(eye(Nfaces),L0),'fro') > tol
%    keyboard
%end

% L0
L0ids = zeros(Nfp,7);
L0vals = zeros(Nfp,7);
if max(sum(abs(L0)>0,2)) > 7
    keyboard
end
for i = 1:Nfp
    tmp = find(L0(i,:));
    L0vals(i,1:length(tmp)) = L0(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 7
        tmp = [tmp zeros(1,7-length(tmp))];
    end
    L0ids(i,:) = tmp;
end

max_EL_nnz = Nfp + 3;
ELids = zeros(Np,max_EL_nnz);
ELvals = zeros(Np,max_EL_nnz);
if max(sum(abs(EL)>0,2)) ~= max_EL_nnz
    keyboard
end
for i = 1:Np
    tmp = find(EL(i,:));
    ELvals(i,1:length(tmp)) = EL(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < max_EL_nnz
        tmp = [tmp zeros(1,max_EL_nnz-length(tmp))];
    end
    ELids(i,:) = tmp;
end

fprintf(fid, '%% BB L0 ids\n');
for n=1:Nfp
    for m=1:7
        fprintf(fid, '%d ', L0ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB L0 values\n');
for n=1:Nfp
    for m=1:7
        fprintf(fid, '%17.15E ', L0vals(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB EL ids\n');
for n=1:Np
    for m=1:max_EL_nnz
        fprintf(fid, '%d ', ELids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB EL values\n');
for n=1:Np
    for m=1:max_EL_nnz
        fprintf(fid, '%17.15E ', ELvals(n,m));
    end
    fprintf(fid, '\n');
end

%degree raise/lower operators along traces
BBRaise = bern_basis_tri(N+1,r2Dp1,s2Dp1)\bern_basis_tri(N,r2Dp1,s2Dp1);
BBRaise(abs(BBRaise)<tol) = 0;

Nfpp1 = (N+2)*(N+3)/2;
BBRaiseIds  = zeros(Nfpp1,3);
BBRaiseVals = zeros(Nfpp1,3);

for i = 1:Nfpp1
    tmp = find(BBRaise(i,:));
    BBRaiseVals(i,1:length(tmp)) = BBRaise(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 3
        tmp = [tmp zeros(1,3-length(tmp))];
    end
    BBRaiseIds(i,:) = tmp;
end

VB2D = bern_basis_tri(N,r2D,s2D);
V2D = Vandermonde2D(N, r2D,s2D);

Nm1 = N-1;
if (N>1)
    [r2Dm1 s2Dm1] = Nodes2D(Nm1); [r2Dm1 s2Dm1] = xytors(r2Dm1,s2Dm1);
    VB2Dm1 = bern_basis_tri(Nm1,r2Dm1,s2Dm1);
    V2Dm1 = Vandermonde2D(Nm1, r2Dm1,s2Dm1);
else
    r2Dm1 =0;
    s2Dm1 =0;
    VB2Dm1 = 1;
    V2Dm1 = 1;
end;

Nfpm1 = (Nm1+1)*(Nm1+2)/2;

BBLower = V2D\VB2D;
BB = [];
sk =1;
for n=0:N
    for m=0:N-n
        if n+m<N
            BB = [BB; BBLower(sk,:)];
        end
        sk = sk+1;
    end
end
BBLower = VB2Dm1\V2Dm1*BB;

fprintf(fid, '%% BB degree raise ids\n');
for n=1:Nfpp1
    for m=1:3
        fprintf(fid, '%d ', BBRaiseIds(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB degree raise values\n');
for n=1:Nfpp1
    for m=1:3
        fprintf(fid, '%17.15E ', BBRaiseVals(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB degree lower matrix\n');
for n=1:Nfpm1
    for m=1:Nfp
        fprintf(fid, '%17.15E ', BBLower(n,m));
    end
    fprintf(fid, '\n');
end


%% elliptic patch problem
K = 5;

VX = [-1, 1,      0,          0,           0,         5/3,        -5/3,         0];
VY = [ 0, 0,sqrt(3),  1/sqrt(3),-7*sqrt(3)/9, 8*sqrt(3)/9, 8*sqrt(3)/9, 1/sqrt(3)];
VZ = [ 0, 0,      0,2*sqrt(6)/3, 4*sqrt(6)/9, 4*sqrt(6)/9, 4*sqrt(6)/9,-2*sqrt(6)/3];

EToV = [1,2,3,4;
    1,2,4,5;
    2,3,4,6;
    3,1,4,7;
    1,3,2,8];

BCType = [0,0,0,0;
    0,0,0,0;
    0,0,0,0;
    0,0,0,0;
    0,0,0,0];

StartUp3D;

% build weak Poisson operator matrices
[A, M] = PoissonIPDG3D();

%% hack since we know W&B face 1 nodes are first
vmapP = reshape(vmapP, Nfp*Nfaces, K);
idsP = vmapP(:,1);
subind = [(1:Np)';idsP];

subA = full(A(subind,subind));
subM = full(M(subind,subind));

spy(abs(subA)>1e-10);
condSubA = cond(subA);

[B,d] = eig(subA, subM);

%% A = S + lambda*M
%%   = M*(M\S + lambda*I)
%%   ~ J*Mhat*(Mhat\Shat/hscale2 + lambda*I)
%%   ~ J*Mhat*Bhat*(d/scale + lambda*I)/Bhat
%% inv(A) ~ Bhat*inv(J*(d/scale+lambda*I))*inv(Mhat*Bhat)

forwardMatrix = inv(subM*B);
diagOp = diag(d);
backwardMatrix = B;

NpP = size(subA,1);
fprintf(fid, '%% stencil size for IPDG OAS NpP\n');
fprintf(fid, '%d\n', NpP);
fprintf(fid, '%% forward matrix [ change of basis for IPDG OAS precon ]\n');
for n=1:NpP
    for m=1:NpP
        fprintf(fid, '%17.15E ', forwardMatrix(n,m));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '%% diagOp \n');
for n=1:NpP
    fprintf(fid, '%17.15E ', diagOp(n));
    fprintf(fid, '\n');
end
fprintf(fid, '%% backward matrix [ reverse change of basis for IPDG OAS precon ]\n');
for n=1:NpP
    for m=1:NpP
        fprintf(fid, '%17.15E ', backwardMatrix(n,m));
    end
    fprintf(fid, '\n');
end

%% permutations
mVXYZ = [-1,+1,-1,-1;
	 -1,-1,+1,-1;
	 -1,-1,-1,+1];
Nverts = 4;
cnt = 0;
for v1=1:Nverts
  for v2=1:Nverts
    for v3=1:Nverts
      for v4=1:Nverts

	if(v1~=v2 & v1~=v3 & v1~=v4 & v2~=v3 & v2 ~=v4 & v3~=v4)
	  vX1 = mVXYZ(:,v1)';
	  vX2 = mVXYZ(:,v2)';
	  vX3 = mVXYZ(:,v3)';
	  vX4 = mVXYZ(:,v4)';
				  
	  permRST = -0.5*(1+r+s+t)*vX1+0.5*(1+r)*vX2+0.5*(1+s)*vX3+0.5*(1+t)*vX4;

	  [v1,v2,v3,v4]
	  
	  permr = permRST(:,1);
	  perms = permRST(:,2);
	  permt = permRST(:,3);

	  for n=1:Np
	    for m=1:Np
	      dist(n,m) = (r(n)-permr(m))^2 + (s(n)-perms(m))^2 + (t(n)-permt(m))^2;
	    end
	  end
	  [foo,ids] = min(dist);
	  cnt = cnt+1;
	  ids
	  pmap(cnt, :) = [v1,v2,v3,v4,ids];

	end
      end
    end
  end
end
pmap

%% degree raising interpolation
[rP1,sP1,tP1] = Nodes3D(N+1);
[rP1,sP1,tP1] = xyztorst(rP1,sP1,tP1);

VP1 = Vandermonde3D(N, rP1, sP1, tP1);
IP1 = VP1/V;
NpP1 = length(rP1);

fprintf(fid, '%% degree raising interpolation matrix\n')
fprintf(fid, '%d %d\n', NpP1, Np);
for n=1:NpP1
  for m=1:Np
    fprintf(fid, '%17.15E ', IP1(n,m));
  end
  fprintf(fid, '\n');
end

%% degree lowering interpolation
if(N>1)
  [rM1,sM1,tM1] = Nodes3D(N-1);
[rM1,sM1,tM1] = xyztorst(rM1,sM1,tM1);
else
%% hard code degree 0
rM1 = -1/2; % -1-1-1+1/4
sM1 = -1/2;
tM1 = -1/2;
end

VM1 = Vandermonde3D(N, rM1, sM1, tM1);
IM1 = VM1/V;
NpM1 = length(rM1);

fprintf(fid, '%% degree lowering interpolation matrix\n')
fprintf(fid, '%d %d\n', NpM1, Np);
for n=1:NpM1
  for m=1:Np
    fprintf(fid, '%17.15E ', IM1(n,m));
  end
  fprintf(fid, '\n');
end




fclose(fid);
