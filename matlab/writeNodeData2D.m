
function writeNodeData2D(inN)

Globals2D;

N = inN;
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);

Np = length(r);
Nfp = N+1;
Nfaces = 3;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)';
faceNodes2   = find( abs(r+s) < NODETOL)';
faceNodes3   = find( abs(r+1) < NODETOL)';
FaceNodes  = [faceNodes1;faceNodes2;faceNodes3]';

V = Vandermonde2D(N, r, s);
[Dr,Ds] = Dmatrices2D(N, r, s, V);
LIFT = Lift2D(N, FaceNodes, r, s);

fname = sprintf('triangleN%02d.dat', N);

fid = fopen(fname, 'w');

fprintf(fid, '%% degree N\n');
fprintf(fid, '%d\n', N);
fprintf(fid, '%% number of nodes\n');
fprintf(fid, '%d\n', Np);
fprintf(fid, '%% node coordinates\n');
for n=1:Np
    fprintf(fid, '%17.15E %17.15E\n', r(n), s(n));
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

MM = inv(transpose(V))/V;
fprintf(fid, '%% reference mass matrix\n');
for n=1:Np
    for m=1:Np
        fprintf(fid, '%17.15E ', MM(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% FaceNodes\n');
for f=1:Nfaces
    for m=1:Nfp
        fprintf(fid, '%d ', FaceNodes(m,f)-1); %% adjust for 0-indexing
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
[plotR,plotS] = EquiNodes2D(N+4);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunay(plotR,plotS)-1;

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1);

%% transform to bi-unit triangle
[plotR,plotS] = xytors(plotR,plotS);

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = Vandermonde2D(N, plotR,plotS)/V;

%% output plot nodes
fprintf(fid, '%% number of plot nodes\n');
fprintf(fid, '%d\n', plotNp);
fprintf(fid, '%% plot node coordinates\n');
for n=1:plotNp
    fprintf(fid, '%17.15E %17.15E\n', plotR(n), plotS(n));
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
    fprintf(fid, '%d %d %d\n' ,...
        plotEToV(n,1),plotEToV(n,2),plotEToV(n,3));
end

%% volume cubature
[cubr,cubs,cubw] = Cubature2D(3*N);
cInterp = Vandermonde2D(N, cubr, cubs)/V;
Np
Ncub = length(cubr)

fprintf(fid, '%% number of volume cubature nodes\n');
fprintf(fid, '%d\n', length(cubr));

fprintf(fid, '%% cubature node coordinates and weights\n');
for n=1:Ncub
	fprintf(fid, '%17.15E %17.15E %17.15E\n', cubr(n), cubs(n), cubw(n));
end

fprintf(fid, '%% cubature interpolation matrix\n');
for n=1:Ncub
    for m=1:Np
        fprintf(fid, '%17.15E ', cInterp(n,m));
    end
    fprintf(fid, '\n');
end

cV = Vandermonde2D(N, cubr, cubs);
cV'*diag(cubw)*cV;

[cVr,cVs] = GradVandermonde2D(N, cubr, cubs);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)
if 0
  Ncut = 0;   scut = 2; 
 filter = Filter2D(V,N,Ncut,scut);    
 cubFilterProject  = filter*cubProject;  
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

fprintf(fid, '%% cubature projection matrix\n');
for n=1:Np
    for m=1:Ncub
        fprintf(fid, '%17.15E ', cubProject(n,m));
    end
    fprintf(fid, '\n');
end
cubProject
testIdentity = cubProject*cInterp;

if 0
fprintf(fid, '%% cubature filter projection matrix\n');
for n=1:Np
    for m=1:Ncub
        fprintf(fid, '%17.15E ', cubFilterProject(n,m));
    end
    fprintf(fid, '\n');
end
end
%% surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
%z = JacobiGL(0,0,N);
%zV = Vandermonde1D(N,z);
%w = sum(inv(zV*transpose(zV)));

Nfi = length(z);

ir = [z,-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, -z];
iw = [w,w,w];

sV = Vandermonde2D(N, ir(:), is(:));
sInterp = sV/V;

iInterp = [sInterp(1:Nfi,FaceNodes(:,1));sInterp(Nfi+1:2*Nfi,FaceNodes(:,2));sInterp(2*Nfi+1:3*Nfi,FaceNodes(:,3))];

fprintf(fid, '%% number of surface integration nodes per face\n');
fprintf(fid, '%d\n', length(z));
fprintf(fid, '%% surface integration interpolation matrix\n');
for n=1:Nfi*Nfaces
    for m=1:Nfp
        fprintf(fid, '%17.15E ', iInterp(n,m));
    end
    fprintf(fid, '\n');
end

bInterp = [];
bInterp(1:Nfi,1:Nfp) = iInterp(1:Nfi,:);
bInterp(Nfi+1:2*Nfi,Nfp+1:2*Nfp) = iInterp(Nfi+1:2*Nfi,:);
bInterp(2*Nfi+1:3*Nfi,2*Nfp+1:3*Nfp) = iInterp(2*Nfi+1:3*Nfi,:);

%% integration node lift matrix
iLIFT = V*V'*sInterp'*diag(iw(:));
size(iLIFT)
size(iInterp)
max(max(abs(iLIFT*bInterp-LIFT)))

fprintf(fid, '%% surface integration lift matrix\n');
for n=1:Np
    for m=1:Nfi*Nfaces
        fprintf(fid, '%17.15E ', iLIFT(n,m));
    end
    fprintf(fid, '\n');
end

cubDrT*ones(Ncub,1)
nr = [zeros(Nfi,1);ones(Nfi,1);-ones(Nfi,1)];
ns = [-ones(Nfi,1);ones(Nfi,1);zeros(Nfi,1)];
sJ = [ones(Nfi,1);ones(Nfi,1);ones(Nfi,1)];
cubDrT*ones(Ncub,1) - iLIFT*(nr.*sJ)
cubDsT*ones(Ncub,1) - iLIFT*(ns.*sJ)


%% BB

addpath('./bern')
tol=1e-6;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
Np = (N+1)*(N+2)/2;
Nfaces = 3;
Nfp = N+1;

r1D = JacobiGL(0,0,N);
[r1Dq w1Dq] = JacobiGQ(0,0,N);
[rq sq wq] = Cubature2D(3*N);
V1Dq = bern_basis_1D(N,r1Dq);

V = Vandermonde2D(N,r,s);
[VB Vr Vs V1 V2 V3] = bern_basis_tri(N,r,s);
VB1D = bern_basis_1D(N,r1D);
invVB = inv(VB);

%% write VDM for conversion

% VB(abs(VB)<tol) = 0;
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

D1 = VB\V1; D1(abs(D1)<tol) = 0;
D2 = VB\V2; D2(abs(D2)<tol) = 0;
D3 = VB\V3; D3(abs(D3)<tol) = 0;

D1ids = zeros(Np,3);
D2ids = zeros(Np,3);
D3ids = zeros(Np,3);
Dvals = zeros(Np,3);
D2vals = zeros(Np,3);
D3vals = zeros(Np,3);
for i = 1:Np
    tmp = find(D1(i,:));
    Dvals(i,1:length(tmp)) = D1(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 3
        tmp = [tmp zeros(1,3-length(tmp))];
    end
    D1ids(i,:) = tmp; 
   
    
    tmp = find(D2(i,:));
    D2vals(i,1:length(tmp)) = D2(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 3
        tmp = [tmp zeros(1,3-length(tmp))];
    end
    D2ids(i,:) = tmp;
    
    tmp = find(D3(i,:));
    D3vals(i,1:length(tmp)) = D3(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 3
        tmp = [tmp zeros(1,3-length(tmp))];
    end
    D3ids(i,:) = tmp;
end

fprintf(fid, '%% sparse BB D1 matrix ids\n');
for n=1:Np
    for m=1:3
        fprintf(fid, '%d ', D1ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D2 matrix ids\n');
for n=1:Np
    for m=1:3
        fprintf(fid, '%d ', D2ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D3 matrix ids\n');
for n=1:Np
    for m=1:3
        fprintf(fid, '%d ', D3ids(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% sparse BB D1,D2,D3 matrix values\n');
for n=1:Np
    for m=1:3
        fprintf(fid, '%17.15E ', Dvals(n,m));
    end
    fprintf(fid, '\n');
end

%% cubature

VBq = bern_basis_tri(N,rq,sq);
M = VBq'*diag(wq)*VBq;
PBq = M\(VBq'*diag(wq));

if (norm(PBq*VBq-eye(Np),'fro')>tol)
    keyboard
end
fprintf(fid, '%% BB cubature interpolation matrix\n');
for n=1:Ncub
    for m=1:Np
        fprintf(fid, '%17.15E ', VBq(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB cubature projection matrix\n');
for n=1:Np
    for m=1:Ncub
        fprintf(fid, '%17.15E ', PBq(n,m));
    end
    fprintf(fid, '\n');
end

%% LIFT

fmask1   = find( abs(s+1) < tol)';
fmask2   = find( abs(r+s) < tol)';
fmask3   = find( abs(r+1) < tol)';
Fmask  = [fmask1;fmask2;fmask3]';

Emat = zeros(Np, Nfaces*Nfp);
% face 1
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(N, faceR);
massEdge1 = inv(V1D*V1D');
Emat(Fmask(:,1),1:Nfp) = massEdge1;
% face 2
faceR = r(Fmask(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massEdge2;
% face 3
faceS = s(Fmask(:,3));
V1D = Vandermonde1D(N, faceS);
massEdge3 = inv(V1D*V1D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;
% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
% convert to BB
LIFT = VB\(LIFT * kron(eye(Nfaces),VB1D));

EL1 = [];
for j = 0:N
    l(j+1) = (-1)^j * nchoosek(N,j)/(1+j);
    Ej = bern_basis_1D(N,r1Dq)\bern_basis_1D(N-j,r1Dq);
    EL1 = [EL1; l(j+1)*Ej'];
end
EL1(abs(EL1)<tol) = 0;

ENp1 = bern_basis_1D(N+1,JacobiGL(0,0,N+1))\bern_basis_1D(N,JacobiGL(0,0,N+1));

L0 = (N+1)^2/2 * ENp1'*ENp1;

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
if norm(LIFT-EL*kron(eye(Nfaces),L0),'fro') > tol
    keyboard
end

% L0 vals
L0vals = zeros(Nfp,3);
for i = 1:Nfp
    if (i==1)
        L0vals(i,2:3) = [L0(i,i) L0(i,i+1)]; %tridiagonal fix
    elseif (i==Nfp)
        L0vals(i,1:2) = [L0(i,i-1) L0(i,i)];
    else
        L0vals(i,:) = [L0(i,i-1) L0(i,i) L0(i,i+1)];
    end
end

max_EL_nnz = Nfp + 2;
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


fprintf(fid, '%% BB L0 values\n');
for n=1:Nfp
    for m=1:3
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
BBRaise = bern_basis_1D(N,r1Dq)\bern_basis_1D(N-1,r1Dq);
BBRaise(abs(BBRaise)<tol) = 0;

BBRaiseIds  = zeros(N+1,2);
BBRaiseVals = zeros(N+1,2);

for i = 1:N+1
    tmp = find(BBRaise(i,:));
    BBRaiseVals(i,1:length(tmp)) = BBRaise(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 2
        tmp = [tmp zeros(1,2-length(tmp))];
    end
    BBRaiseIds(i,:) = tmp;
end

[r1Dp1] = JacobiGQ(0,0,N+1);
VB1Dp1 = bern_basis_1D(N+1,r1Dp1);
V1Dp1 = Vandermonde1D(N+1, r1Dp1);

[r1D] = JacobiGQ(0,0,N);
VB1D = bern_basis_1D(N,r1D);
V1D = Vandermonde1D(N, r1D);

BBLower = V1Dp1\VB1Dp1;
BBLower = VB1D\V1D*BBLower(1:N+1,1:N+2);

fprintf(fid, '%% BB degree raise ids\n');
for n=1:N+1
    for m=1:2
        fprintf(fid, '%d ', BBRaiseIds(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB degree raise values\n');
for n=1:N+1
    for m=1:2
        fprintf(fid, '%17.15E ', BBRaiseVals(n,m));
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%% BB degree lower matrix\n');
for n=1:N+1
    for m=1:N+2
        fprintf(fid, '%17.15E ', BBLower(n,m));
    end
    fprintf(fid, '\n');
end



%% elliptic patch problem
K = 4;

%VX = [-1,+1,-1,-1,+1,-3];
%VY = [-1,-1,+1,-3,+1,-1];

VX = [-1,1,0,0,2,-2];
VY = [0,0,sqrt(3),-sqrt(3),sqrt(3),sqrt(3)];

EToV = [1,2,3;
	2,1,4;
	3,2,5;
	1,3,6];

BCType = [0,0,0;
	  0,0,0;
	  0,0,0;
	  0,0,0];

StartUp2D;


% choose order to integrate exactly
Nint = ceil(2*N/2);

% build cubature nodes for all elements
CubatureOrder = 2*(Nint+1); 
cub = CubatureVolumeMesh2D(CubatureOrder);
  
% build Gauss node data for all element faces
NGauss = (Nint+1); 
gauss = GaussFaceMesh2D(NGauss);
  
% build weak Poisson operator matrices
[A, M] = CurvedPoissonIPDG2D();

%% hack since we know W&B face 1 nodes are first
vmapP = reshape(vmapP, Nfp*Nfaces, K);
idsP = vmapP(:,1)
subind = [(1:Np)';idsP];

subA = full(A(subind,subind))
subM = full(M(subind,subind));

spy(abs(subA)>1e-10);
condSubA = cond(subA);

[B,d] = eig(subA, subM)

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



fclose(fid)
