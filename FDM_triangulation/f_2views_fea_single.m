function [ X ] = f_2views_fea_single( PArray, Img, c, d )

    nCam = length(d);
    
    T = 10;
    t = 0;
    while t<T
        t = t+1;
        minsub = randsample(nCam,2);
        P1 = PArray(3*(minsub(1)-1)+1:3*(minsub(1)-1)+3,:); 
        P2 = PArray(3*(minsub(2)-1)+1:3*(minsub(2)-1)+3,:); 
        u1 = [Img(2*(minsub(1)-1)+1:2*(minsub(1)-1)+2) ; 1]; 
        u2 = [Img(2*(minsub(2)-1)+1:2*(minsub(2)-1)+2) ; 1] ; 
        feaX = intsec2views_midpoint(P1,P2,u1,u2);
        feaX = feaX(1:3);
        fea = c * feaX + d > 0;
        
        if nnz(fea) == nCam
            X = feaX;
            break;
        else
            continue;
        end
    end
    
    if ~exist('X','var')
        X = f_init( PArray, Img); 
    end
    
end



