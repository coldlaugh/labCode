function [B,B_broken,G,relaxation]=break_bonds(X,X0,B,B_broken,S,G,W)
    
    BC = @(x) [x(1)-round(x(1)/W)*W,x(2)];
    r = zeros(size(B,1),1);
    for i = 1:size(B)
        n1 = B(i,1); n2 = B(i,2);
        if any(S==n1) || any(S==n2)
            r(i) = norm(BC(X(n1,:) - X(n2,:))) - norm(BC(X0(n1,:) - X0(n2,:)));
        else
            r(i) = 0;
        end
    end

    [maxforce maxforcebond] = max(r./G);
    relaxation = false;
    if (maxforce > 1)
        B_broken = [B_broken;
            B(maxforcebond,:)];
        B(maxforcebond,:) = [];
        G(maxforcebond)=[];
        relaxation = true;
    end  
end