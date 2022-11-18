function [ XArray ] = f_2views_fea( PArray, ImgArray )

N = size(ImgArray,2);
M = size(PArray,1)/3;

XArray = zeros(3,N);
[A1,~,~,~,C,D] = f_gen_coef_2(PArray, ImgArray);

for n = 1:N

    a1 = A1(M*(n-1)+1:M*n,:);
    c  = C;
    d  = D;
    
    fid = find(isfinite(sum(a1,2)));
    c  = c (fid,:);
    d  = d (fid,:);
    
    id1 = (fid-1)*3+1;
    id2 = (fid-1)*3+2;
    id3 = (fid-1)*3+3;
    pid = [id1';id2';id3'];
    pid = pid(:);
    id1 = (fid-1)*2+1;
    id2 = (fid-1)*2+2;
    iid = [id1';id2'];
    iid = iid(:);
    
    XArray(:,n) = f_2views_fea_single( PArray(pid,:), ImgArray(iid,n), c, d );  
end % end of for

end



