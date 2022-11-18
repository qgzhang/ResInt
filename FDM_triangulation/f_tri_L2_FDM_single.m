function [ X, mres, nitr] = f_tri_L2_FDM_single(a1, a2, b1, b2, c, d, X0, aeps)


    mres0 = Inf;
    mres = 1e10;
    flag = true;

    tol = 1e-8;  % 1e-16
    tol = 1e-4;
    maxitr = 100; % 100
    nitr = 0;
    X = X0;
    
    while abs(mres0-mres)> tol  && nitr<maxitr
        nitr = nitr+1;
        % find the active constraints

%         res = ( (a1*X+b1).*(a1*X+b1) + (a2*X+b2).*(a2*X+b2) ) ./ ( (c*X+d).*(c*X+d) );
        res = sqrt( (a1*X+b1).^2 + (a2*X+b2).^2 ) ./ (c*X+d);
        %tt = [res(1:M),res(M+1:2*M),res(2*M+1:3*M),res(3*M+1:4*M)];
        mres0 = mres;
        mres = max(res);
        if mres>mres0
            flag = false;
            break;
        end
        Active = find(mres-res<=mres*aeps)';
%         if length(Active)>4
%             fprintf('No. Active is: %d\n', length(Active))
%         end
        
        [grad] = f_grad_L2(a1(Active,:), a2(Active,:), b1(Active), b2(Active), c(Active,:), d(Active), X);
    %    grad = sparse(grad);
    %     fprintf('improvement = %.10f\n', );
        % get a descent direction
        % opt: 1 - max dot product; 2 - feasibility test
        

        %[ r, dir ] = f_finddir( -grad); 
        [ dir, ~ ] = f_meb( -grad);  r = norm(dir);
        
%         [ r, dir ] = f_finddir_poly( -grad); 
    
%         if first==1 % for the guage freedom: fix X1(3) = 1
%             dir(3) = 0;
%         end
    %     if mod(nitr,10)==0
%             fprintf('%d: nnz(dir)-%d, impv-%.10f\n', nitr, nnz(dir), mres_1000-mres);
%             mres_1000 = mres;
    %     end
        if r<tol % none valid descent direction
            break;
        end
        dir = dir/norm(dir);
        
%         fprintf('(%d) dir = %.6f %.6f %.6f\n', nitr, dir(1), dir(2), dir(3));
        % get a step size
        % opt: 1 - incremental; 2 - bisection
        [lb ub] = f_getbound(c,d,X,dir);
        [ alpha ] = f_stepsize_L2( a1, a2, b1, b2, c, d, X, dir, Active, mres, lb, ub, 1 );
%         foo(a1, a2, b1, b2, c, d, X, dir, Active, alpha, max(lb), min(1.1*alpha,ub));
%         [miny,mina] = foo2(a1, a2, b1, b2, c, d, X, dir, Active, alpha)
%         clf
        if alpha < tol
            break;
        end
        % update X
        X0 = X;
        X = X + alpha * dir;

    end % while

    if ~flag
        X = X0;
        mres = mres0;
    end

%     res_cons = (aa*X+bb)./(cc*X+dd);
%     mres = max(res_cons); % max error (converged error)


end


