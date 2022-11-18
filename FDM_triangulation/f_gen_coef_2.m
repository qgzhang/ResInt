function [A1,A2,B1,B2,C,D] = f_gen_coef_2(PArray, ImgArray)

    M = size(PArray,1)/3;
    N = size(ImgArray,2);

%     A1 = zeros(M*N,3);
%     A2 = zeros(M*N,3);
%     B1 = zeros(M*N,1);
%     B2 = zeros(M*N,1);
%     C  = zeros(M,3);
%     D  = zeros(M,1);
    
    R1 = PArray(1:3:end,1:3);
    R2 = PArray(2:3:end,1:3);
    R3 = PArray(3:3:end,1:3);
    T1 = PArray(1:3:end,4);
    T2 = PArray(2:3:end,4);
    T3 = PArray(3:3:end,4);
    
    U = ImgArray(1:2:end,:);
    V = ImgArray(2:2:end,:);
    
    A1 = repmat(R1,N,1) - repmat(U(:),1,3).*repmat(R3,N,1);
    A2 = repmat(R2,N,1) - repmat(V(:),1,3).*repmat(R3,N,1);
    B1 = repmat(T1,N,1) - U(:).*repmat(T3,N,1);
    B2 = repmat(T2,N,1) - V(:).*repmat(T3,N,1);
    C = R3;
    D = T3;
    
    
end % end of function gen_coef()

