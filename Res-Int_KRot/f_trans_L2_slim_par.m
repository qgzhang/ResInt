function [ TArray, mres, nitr] = f_trans_L2_slim_par(A, u, XArray, initT, aeps) %done

N = size(u{1},2);
M = length(A);
TArray = zeros(3,M);
mres = zeros(1,M);

nitr = 0;

if exist('initT','var') && ~isempty(initT)
    T0 = initT;
else
    T0 = ones(3,M);
end

parfor m = 1:M
    
    R = A{m};
    Img = u{m};
    
    fid = isfinite(Img(1,:));
    X = XArray(:,fid);
    Img = Img(:,fid);
    n = size(X,2);
    
    R1 = R(1,:);
    R2 = R(2,:);
    R3 = R(3,:);
    X1 = X(1,:);
    X2 = X(2,:);
    X3 = X(3,:);
    U = Img(1,:);
    V = Img(2,:);
    b1 = ( repmat(R1,n,1) - repmat(U(:),1,3).*repmat(R3,n,1) ) ;
    B11 = b1(:,1);
    B12 = b1(:,2);
    B13 = b1(:,3);
    b1 = B11.*X1(:) + B12.*X2(:) + B13.*X3(:);
    b2 = ( repmat(R2,n,1) - repmat(V(:),1,3).*repmat(R3,n,1) ) ;
    B21 = b2(:,1);
    B22 = b2(:,2);
    B23 = b2(:,3);
    b2 = B21.*X1(:) + B22.*X2(:) + B23.*X3(:);
    d = repmat(R3,n,1);
    D1 = d(:,1);
    D2 = d(:,2);
    D3 = d(:,3);
    d = D1.*X1(:) + D2.*X2(:) + D3.*X3(:);
   
    a1 = [ones(n,1), zeros(n,1), -U(:)];
    a2 = [zeros(n,1), ones(n,1), -V(:)];
    c  = repmat([0,0,1], n , 1);

    % init T
    Tm0 = T0(:,m);
     % number of violation
    if Tm0(3) <= max(-d) + 1e-10
        Tm0(3) = max(-d) + 1;
    end
  
    [T, mres_s, nitr_s] = c_tri_single_L2(a1, a2, b1, b2, c, d, Tm0, aeps);


    nitr = nitr + nitr_s;
    TArray(:,m) = T;
    mres(m) = mres_s;

end % for n=1:N

end



