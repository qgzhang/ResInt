function [grad] = f_grad_L2(A1, A2, B1, B2, C, D, X) % done
% grad: [grad1, grad2, ... ];
% A1, ..., D: active a1 to d

% quotient rule: (f/g)' = (f'g - fg')/g^2
% since only the direction matters, we can neglect the positive denominator g^2


[active, d] = size(A1);
grad = zeros(d, active);
GX = (A1*X+B1).^2 + (A2*X+B2).^2;
for i = 1:active
    
%     grad(:,i) = AA(i,:)' * (CC(i,:)*X + DD(i)) - CC(i,:)' * (AA(i,:)*X+BB(i));
    a1 = A1(i,:);
    a2 = A2(i,:);
    b1 = B1(i,:);
    b2 = B2(i,:);
    c = C(i,:);
    d = D(i,:);
    
    gx = GX(i);
    
    grad(:,i) = (c*X+d)/sqrt(gx) * ((a1*X+b1)*a1' + (a2*X+b2)*a2')  -  sqrt(gx)*c';
    
    grad(:,i) = grad(:,i)/norm(grad(:,i));
    %grad(:,act) = grad(:,act)/ (grad(:,act)'*grad(:,act));
end

% restrict guage freedom: by the paper, X_1(3) = 1, T1 = [0;0;1];
% which means: grad(X(3)) = 0, grad( X(3*N+1: 3*N+3) ) = [0;0;0];
% grad(3,:) = 0;
% grad(3*N+1:3*N+3,:) = repmat([0;0;0],1,size(grad,2));

end
