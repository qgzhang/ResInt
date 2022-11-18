function [lb,ub] = f_getbound(cc, dd, X,dir)

    l_id = find(cc*dir>0); 
    if isempty(l_id)
        lb = 0;
    else
        %lb = max( (-cc(l_id,:)*X) ./ (cc(l_id,:)*dir) );
        lb = max( (-dd(l_id)-cc(l_id,:)*X) ./ (cc(l_id,:)*dir) );
        lb = max(0,lb);
    end

    u_id = find(cc*dir<0);
    if isempty(u_id)
        ub = Inf;
    else
        %ub = min( (-cc(u_id,:)*X) ./ (cc(u_id,:)*dir) );    
        ub = min( (-dd(u_id)-cc(u_id,:)*X) ./ (cc(u_id,:)*dir) );
    end
end

