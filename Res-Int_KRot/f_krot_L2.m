function [X, T, mres, nitr] = f_krot_L2(R, Imgv, A, u, tol, X0, T0,  aeps, parallel, thres)  %done
% X: [X1; X2; ...; Xn; t1; t2; ...; tm]
% aa,cc: reprojection error = aa*X ./ bb*X
mres = 100;
nitr = 0;

N = size(u{1},2);
M = length(A);

if exist('X0','var') && ~isempty(X0)
    X = X0;
else
    X = ones(3,N);
    X0 = X;
end

if exist('T0','var') && ~isempty(T0)
    T = T0;
else
    T = ones(3,M);
    T0 = T;
end

P = [R , T(:)];
hit = 0;
fprintf(' itr: | elapsed_time\n');
while 1
    nitr = nitr+1;
    
    if parallel==0
        [ X, mres_tri, ~] = f_tri_L2_slim(P, Imgv, X0, aeps);
        [ T, mres_trans, ~] = f_trans_L2_slim(A, u, X, T0, aeps);
    else
        [ X, mres_tri, ~] = f_tri_L2_slim_par(P, Imgv, X0, aeps);
        [ T, mres_trans, ~] = f_trans_L2_slim_par(A, u, X, T0, aeps);
    end

    P(:,4) = T(:);

    mres0 = mres;
%         res = f_res(X,P,Img,Norm);
%        mres = max(max(res));
    mres = max(max(mres_trans));

    if mres0-mres < mres*1e-3 || mres0-mres < tol
        hit = hit+1;
        fprintf('%4d: | %8.2f\n', nitr, toc);
    else
        fprintf('%4d: | %8.2f\n', nitr, toc);
        hit = hit - 1;
        hit = max(0,hit);
    end

    if hit>=thres
        %fprintf('(%d) -  %.12f -> %.12f = %.12f\n', nitr, mres0, mres, mres0-mres);
        break;
    end
    
    
    
    T0 = T;
    X0 = X;
end
    
end