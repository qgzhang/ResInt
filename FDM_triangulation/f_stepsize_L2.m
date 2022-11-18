function [ alpha, val, nitr ] = f_stepsize_L2( a1, a2, b1, b2, c, d, X, dir, Active, oldval, lb0, ub0, opt ) % done

% opt: 
% 1 - 'bisection on alpha'
% 2 - 'bisection on r'
% 3 - 'incremental bisection on r'
% 4 - first intersection (Donne's method)

tol1 = 1e-6;
tol2 = 1e-8;
n = size(a1,1);
% get feasible region [lb, lb]
% [lb0,ub0] = f_getbound(c,d,X,dir);
if lb0>=ub0
    alpha = 0;
    val = oldval;
    return;
end
    
nitr = 0;    


% bisection on alpha
if opt==1
    ub = min(ub0,1e1);
    lb = max(0,lb0);
    r = (ub+lb)/2;
    alpha = lb;
    val = oldval;
    while ub-lb > tol1
        r = (ub+lb)/2;
        Xn = X + (r-tol2)*dir; % Xn is on the negative side of Xc
        Xp = X + (r+tol2)*dir;
        Xc = X + (r)*dir;
        vn = max( sqrt( ( (a1*Xn+b1).*(a1*Xn+b1) + (a2*Xn+b2).*(a2*Xn+b2) ) ./ ( (c*Xn+d).*(c*Xn+d) ) ) );
        vp = max( sqrt( ( (a1*Xp+b1).*(a1*Xp+b1) + (a2*Xp+b2).*(a2*Xp+b2) ) ./ ( (c*Xp+d).*(c*Xp+d) ) ) );
        vc = max( sqrt( ( (a1*Xc+b1).*(a1*Xc+b1) + (a2*Xc+b2).*(a2*Xc+b2) ) ./ ( (c*Xc+d).*(c*Xc+d) ) ) );
%         [vn; vc; vp]
        if vn > vc && vc > vp % increase r
            lb = r;
        elseif vn < vc && vc < vp % decrease r
            ub = r;
        else
            alpha = r;
            val = vc;
            break;
        end
    end
    alpha = r;
    val = vc;
end
% 
% % bisection on r
% if opt==2
%     [alpha, val] = linfrac(a1, a2, b1, b2, c, d, X, dir, lb0, ub0, oldval, tol2);
% end
% 
% 
% % incremental bisection on r
% if opt==3
%     cs = randsample(n,2);
%     val0 = Inf;
%     val = 100000;
%     alpha = 0;
% 
%     while abs(val-val0) > tol1
%         nitr = nitr+1;
%         val0 = val;
%         csid = [cs; n+cs; 2*n+cs; 3*n+cs];
%         [alpha, ~] = linfrac(aa(csid,:), bb(csid), cc(csid,:), dd(csid), X, dir, lb0, ub0, oldval, tol2);
%         X1 = X + alpha*dir;
%         [val,q] = max( (aa*X1+bb)./(cc*X1+dd) );
%         q = mod(q,n);
%         if q==0
%             q=n;
%         end
%         if isempty(find(cs==q,1))
%             cs = [cs;q];
%         else
%             break;
%         end
%     end
% end
% 
% 
% 
% % first intersection (Donne's method)
% if opt==4
%     method = 2;
%     a0 = aa(Active,:);
%     b0 = bb(Active,:);
%     c0 = cc(Active,:);
%     d0 = dd(Active,:);
%     act_dev= ( (c0*X+d0).*(a0*dir) - (a0*X+b0).*(c0*dir) ) ./ ((c0*X+d0).^2);
%     % master constraint: the least negative one, so is the maximum one
%     [~,ind] = max(act_dev);
%     master = Active(ind);
%     Roots = f_root(aa,bb,cc,dd,master,1:size(aa,1),X,dir);
%     if method == 1
%         alpha = min(Roots(Roots>1e-10)); 
%     elseif method == 2
%         val = oldval;
%         brt = 0;
%         rts = Roots(Roots>0);
%         for r = 1:length(rts)
%             Xp = X + rts(r)*dir;
%             resp = (aa*Xp+bb)./(cc*Xp+dd);
%             if val > max(resp)
%                 val = max(resp);
%                 brt = rts(r);
%             end
%         end
%         alpha = brt;
%     end
% end
% 
% 
% 
% 
% end % end of function f_stepsize
% 
% 
% 
% % linear fractional programming solver, using bisection
% function [alpha, r] = linfrac(a1, a2, b1, b2, c, d, X, dir, lb0, ub0, oldval, tol)
% % lb0, ub0: feasible region at the beginning oldval
%     upper = oldval;
%     lower = 0;
%     r = oldval;
%     alpha = 0;
%     n = length(d);
%     
%     A0 = (a1*dir).^2 + (a2*dir).^2;
%     B0 = 2*(a1*dir).*(a1*X+b1) + 2*(a2*dir).*(a2*X+b2);
%     C0 = (a1*X+b1).^2 + (a2*X+b2).^2;
%     cdir = (c*dir);
%     cdir2 = (c*dir).^2;
%     cxpd = c*X+d;
%     cxpd2 = (c*X+d).^2;
%     
%     
%     maxlb1 = lb0;
%     minub1 = ub0;
%     fea = true;
%     
%     while upper-lower > tol
%         r = (upper+lower)/2;
%         lb1 = lb0;
%         ub1 = ub0;
%         lb2 = lb0;
%         ub2 = ub0;
%         
%         A = A0 - r^2 * cdir2;
%         B = B0 - 2*r^2 * cdir .* cxpd;
%         C = C0 - r^2 * cxpd2;
%        
%         roots = [(-B+sqrt(B.*B-4*A.*C))./(2*A) , (-B-sqrt(B.*B-4*A.*C))./(2*A)];
% %         roots = roots(:);
% %         for i=1:length(roots) % eliminate imaginary roots 
% %             if ~isreal(roots(i))
% %                 roots(i) = -1;
% %             end
% %         end
% %         roots = reshape(roots(:),[],2);
%         
%         %%% A ~= 0
%         % A > 0
%         id = A > 1e-10;
%         rts = roots(id,:);
%         rts = sort(rts,2,'ascend');
%         for i = 1:length(id)
%             if isreal(rts(i,1)) % has real roots
%                 lb1 = max(lb1, min(rts(i,1)));
%                 ub1 = min(ub1, max(rts(i,1)));
%             else
%                 fea = false;
%                 break;
%             end
%         end
%         if fea==false
%             lower = r;
%             continue;
%         end
%         
%         % A < 0
%         id = A < -1e-10; 
%         rts = roots(id,:);
%         rts = sort(rts,2,'ascend');
%         for i = 1:length(id)
%             if isreal(rts(i,1)) % has real roots
%                 lb2 = max(lb2, max(rts(i,1)));
%                 ub2 = min(ub2, min(rts(i,1)));
%             end
%         end
%         
%         %%% A == 0, B ~= 0
%         % B > 0
%         id = abs(A)<=1e-10 && B > 1e-10;
%         rts = -C(id)/B(id);
%         ub1 = min(ub1, rts);
%         % B < 0
%         id = abs(A)<=1e-10 && B < -1e-10;
%         rts = -C(id)/B(id);
%         lb1 = max(lb1, rts);
%         
%         %%% A == 0, B==0
%         % C>0
%         id = abs(A)<=1e-10 && abs(B)<=1e-10 && C > 1e-10;
%         if ~isempty(id)
%             lower = r;
%             continue;
%         end
%         
%         if (lb2<=ub1 || ub2>=lb1) && (lb<=ub && lb0<=ub0 && lb<=ub0 && ub>=lb0)
%             upper = r;
%             if ub2>=lb1
%                 alpha = ( ub2+lb1 ) /2;
%             elseif lb2<=ub1
%             	alpha = ( lb2+ub1 ) /2;
%             end
%         else
%             lower = r;
%         end
%     end
% 
% end
% 
% % % get feasible region of the step size
% % function [lb,ub] = f_getbound(cc, dd, X,dir)
% % 
% %     l_id = find(cc*dir>0); 
% %     if isempty(l_id)
% %         lb = 0;
% %     else
% %         %lb = max( (-cc(l_id,:)*X) ./ (cc(l_id,:)*dir) );
% %         lb = max( (-dd(l_id)-cc(l_id,:)*X) ./ (cc(l_id,:)*dir) );
% %         lb = max(0,lb);
% %     end
% % 
% %     u_id = find(cc*dir<0);
% %     if isempty(u_id)
% %         ub = Inf;
% %     else
% %         %ub = min( (-cc(u_id,:)*X) ./ (cc(u_id,:)*dir) );    
% %         ub = min( (-dd(u_id)-cc(u_id,:)*X) ./ (cc(u_id,:)*dir) );
% %     end
% % end