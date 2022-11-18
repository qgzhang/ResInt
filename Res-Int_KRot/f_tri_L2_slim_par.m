function [ XArray, mres, nitr, fail] = f_tri_L2_slim_par(PArray, Imgv, initX, aeps) %done

N = length(Imgv);
XArray = zeros(3,N);
mres = zeros(1,N);
M = size(PArray,1)/3;

nitr = 0;
fail = 0;



R1 = PArray(1:3:end,1:3);
R2 = PArray(2:3:end,1:3);
R3 = PArray(3:3:end,1:3);
T1 = PArray(1:3:end,4);
T2 = PArray(2:3:end,4);
T3 = PArray(3:3:end,4);

clear PArray;

parfor n = 1:N
    Img = Imgv{n};
    U = Img(1:2:end);
    V = Img(2:2:end);
    
    a1 = R1 - repmat(U(:),1,3) .* R3;
    a2 = R2 - repmat(V(:),1,3) .* R3;
    b1 = T1 - U .* T3;
    b2 = T2 - V .* T3;
    c = R3;
    d = T3;
    
 
    fid = find(isfinite(sum(a1,2)));
    a1 = a1(fid,:);
    a2 = a2(fid,:);
    c  = c (fid,:);
    b1 = b1(fid,:);
    b2 = b2(fid,:);
    d  = d (fid,:);
   
    Xn = initX(:,n);
    [X, mres_s, nitr_s] = c_tri_single_L2(a1, a2, b1, b2, c, d, Xn, aeps);

    XArray(:,n) = X;
    mres(n) = mres_s;
    nitr = nitr + nitr_s;

end % for n=1:N

end




