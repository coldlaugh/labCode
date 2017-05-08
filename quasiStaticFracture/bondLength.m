function [lRest] = bondLength(X,B,W)  % the length of bonds

    lRest = zeros(size(B,1),1);
    
    for i = 1: size(B,1)
        n1 = B(i,1); n2 = B(i,2);
        if n1>=0 && n2 >= 0
            dr = X(n1,:) - X(n2,:);
            dr(1) = dr(1) - round(dr(1)/W)*W;
            lRest(i) = norm(dr);
        end
    end

end