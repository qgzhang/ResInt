function [c,rad] = f_mycircle(A,B,C)

a = A-C;
b = B-C;

na = norm(a);
nb = norm(b);
aa = na^2*b - nb^2*a;

% rpt = 1000000;

% tic
% for i=1:rpt
nab = norm(cross(a,b));
c_abc = cross(aa, cross(a,b));
% end
% toc

% tic
% for i=1:rpt
% nab2 = (na*na*nb*nb-(a'*b)*(a'*b));
% c_abc2 = aa'*b*a - aa'*a*b;
% end
% toc




rad = na*nb*norm(a-b) / (2*nab);
c = c_abc / (2*nab^2) + C;
end

