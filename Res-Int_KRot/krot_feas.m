function [U,P,dual,s] = krot_feas(u,A,tol,min_depth,max_depth)
% Solves the known rotation feasibility problem
% Olsson, Eriksson, Hartley, Outlier Removal using Duality, CVPR2010
%
% inputs -  u: 1xD cell with image data.
%           u{i} is of size 3xN, where N is the number of observed points.
%           If point j is not observed in image then u{i}(:,j) = NaN.
%
%        -  A: 1xD cell with estimated orientation matrices.
%
%        -  tol: maximal inlier tolerance.
%
%        -  min_depth: minimal allowed depth.
%
%        -  max_depth: maximal allowed depth.
%
% outputs - U: 3xD cell with 3D points
%
%         - P: 1xD cell with camera matrices
%
%         - dual: solution to the dual variables. dual(i) > 0 indicates that U(:,i) might
%           be an outlier if s > 0
%
% (C) 2010 Carl Olsson (calle@maths.lth.se, carl.a.c.olsson@gmail.com)

[a,a0,b,b0,c,c0,Aeq,Beq] = gen_krot(u,A);
[Linfsol,dual,s] = LinfSolverfeas(a,a0,b,b0,c,c0,Aeq,Beq,tol,min_depth,max_depth);
[U,P] = form_str_mot(u,A,Linfsol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,P] = form_str_mot(u,A,sol)
numpts = size(u{1},2);
numcams = length(A);

U = reshape(sol(1:(3*numpts)),[3 numpts]);

tpart = [1;1;1;sol(3*numpts+1:end)];
% P = cell(size(A));
P = zeros(3*size(A,1),4);
for i=1:length(A)
%     P{i} = [A{i} tpart((i-1)*3+1:i*3)];
    P(3*(i-1)+1:3*i,:) = [A{i} tpart((i-1)*3+1:i*3)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,dual,s] = LinfSolverfeas(a,a0,b,b0,c,c0,Aeq,Beq,tol,min_depth,max_depth)

A1 = [-a-tol*c; a-tol*c; -b-tol*c; b-tol*c];
B1 = [a0+tol*c0; -a0+tol*c0; b0+tol*c0; -b0+tol*c0];

A2 = [-c; c];
B2 = [c0-min_depth; max_depth-c0];

A = [A1 -ones(size(A1,1),1);A2 zeros(size(A2,1),1)];
C = [B1;B2];
B = [sparse(size(A1,2),1); 1];
K.l = size(A,1);
pars.eps = 1e-8;
pars.fid = 0;
[X,Y,info] = sedumi(A,-B,C,K,pars);
s = Y(end);
Y = Y(1:end-1);
res = size(a,1);
dual = X(1:res)+X(res+1:2*res)+X(2*res+1:3*res)+X(3*res+1:4*res)+X(4*res+1:5*res)+X(5*res+1:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,a0,b,b0,c,c0,Aeq,Beq] = gen_krot(u,A)
numvar = 3*size(u{1},2)+3*length(A);
numpts = size(u{1},2);
numcams = length(A);

a = [];
b = [];
c = [];
for i = 1:numcams;
    R = A{i};
    p = u{i};
    visible_points = isfinite(p(1,:));
    numres = sum(visible_points);
    
    ptind = find(visible_points');
    pointcoeff = p(1,visible_points)'*R(3,:)-ones(numres,1)*R(1,:);

    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    tcoeff = [-ones(numres,1) zeros(numres,1) p(1,visible_points)'];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newa = sparse(row,col,data,numres,numvar);

    
    ptind = find(visible_points');
    pointcoeff = p(2,visible_points)'*R(3,:)-ones(numres,1)*R(2,:);

    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    tcoeff = [zeros(numres,1) -ones(numres,1) p(2,visible_points)'];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newb = sparse(row,col,data,numres,numvar);

    ptind = find(visible_points');
    pointcoeff = ones(numres,1)*R(3,:);

    pointcol = [(ptind-1)*3+1 (ptind-1)*3+2 ptind*3];
    pointrow = [1:numres]'*[1 1 1];
    
    tcoeff = [zeros(numres,1) zeros(numres,1) ones(numres,1)];
    tcol = ones(numres,1)*[numpts*3+[(i-1)*3+1:i*3]];
    trow = pointrow;
    
    data = [pointcoeff(:); tcoeff(:)];
    row = [pointrow(:);  trow(:)];
    col = [pointcol(:); tcol(:)];
    newc = sparse(row,col,data,numres,numvar);

    a = [a;newa];
    b = [b;newb];
    c = [c;newc];
    
end

a0 = zeros(size(a,1),1);
b0 = zeros(size(b,1),1);
c0 = zeros(size(c,1),1);

a = a(:,[1:3*numpts 3*numpts+4:end]);
b = b(:,[1:3*numpts 3*numpts+4:end]);
c = c(:,[1:3*numpts 3*numpts+4:end]);

Aeq = zeros(3,numpts*3);
for i = 1:length(A)
    Aeq = [Aeq -inv(A{i})];
end
Beq = zeros(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%