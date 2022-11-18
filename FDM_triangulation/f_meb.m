function [ c, rad ] = f_meb( P )

% return a vector 'dir' that has positive dot product with any pi in P

% P = [p1, p2, ..., pn],  p1 = [p11; p12; ...; p1d]
% dir: new direction
% r: maximum dot product
% opt: 1 - max dotproduct, 2 - feasibility test

% if isempty(P)
%     error('grad is empty\n');
% end

n = size(P,2);



if n==1
    c = P;
    rad = norm(P);
elseif n==2
    c = (P(:,1)+P(:,2))/2;
    rad = norm(c-P(:,1));
elseif n==3
    f = 0;
    seq = [1 2 3; 
           1 3 2; 
           2 3 1];
    for i = 1:n
        pu = P(:,seq(i,1));
        pv = P(:,seq(i,2));
        pr = P(:,seq(i,3));
        c = (pu+pv)/2;
        rad = norm(c-pu);
        if norm(c-pr) <= rad
            f = 1;
            break;
        end
    end
    if f==0
%         [cs,rad,~,~] = circlefit3d(pu',pv',pr'); c = cs';
        [c,rad] = f_mycircle(pu,pv,pr);
    end
    
else
    f = 0;
    n = 4;
    P = P(:,1:n);
    seq = [ 1     2     3 4;
            1     2     4 3;
            1     3     4 2;
            2     3     4 1];
    for i = 1:n
        pu = P(:,seq(i,1));
        pv = P(:,seq(i,2));
        pw = P(:,seq(i,3));
        pr = P(:,seq(i,4));
        [c,rad] = f_meb([pu,pv,pw]);
%         sqrt(sum(([c,c,c]-[pu,pv,pw]).^2))
%         norm(c-pr)
        if norm(c-pr) <= rad
            f = 1;
            break;
        end
    end
    if f==0
        c = [0;0;0];
        rad = 1;
    end
end

end






