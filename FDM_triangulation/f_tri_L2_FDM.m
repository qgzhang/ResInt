function [ XArray, mres, nitr, fail] = f_tri_L2_FDM(PArray, ImgArray, initX, aeps) %done

N = size(ImgArray,2);
XArray = zeros(3,N);
mres = zeros(1,N);
M = size(PArray,1)/3;

nitr = 0;
fail = 0;
[A1,A2,B1,B2,C,D] = f_gen_coef_tri(PArray, ImgArray);

for n = 1:N
    a1 = A1(M*(n-1)+1:M*n,:);
    a2 = A2(M*(n-1)+1:M*n,:);
    b1 = B1(M*(n-1)+1:M*n,:);
    b2 = B2(M*(n-1)+1:M*n,:);

    c  = C;
    d  = D;
    
    fid = find(isfinite(sum(a1,2)));
    a1 = a1(fid,:);
    a2 = a2(fid,:);
    c  = c (fid,:);
    b1 = b1(fid,:);
    b2 = b2(fid,:);
    d  = d (fid,:);
    

    Xn = initX(:,n);

    %fprintf('Scene point %d//%d\n', n, N)
    [X, mres_s, nitr_s] = f_tri_L2_FDM_single(a1, a2, b1, b2, c, d, Xn, aeps);
    %fprintf('\n')
    XArray(:,n) = X;
    mres(n) = mres_s;
    nitr = nitr + nitr_s;

end % for n=1:N

end









