function [B,B_broken,G,relaxation]=break_bonds(X,lRest,B,B_broken,S,G,W)
    
    
    lCurrent = bondLength(X,B,W);
    r = zeros(size(B,1),1);
    for i = 1:size(B)
        n1 = B(i,1); n2 = B(i,2);
        if any(S==n1) || any(S==n2)
            r(i) = lCurrent(i) - lRest(i);
        else
            r(i) = 0;
        end
    end

%     [maxforce maxforcebond] = max(r./G);
    forceRatio = r ./ G ./rRest; % the ratio between the relative extension and the threshold
    relaxation = any(forceRatio >= 1);
    if (maxforce > 1)
        B_broken = [B_broken;
            B(maxforcebond,:)];
        B(maxforcebond,:) = [];
        G(maxforcebond)=[];
        relaxation = true;
    end  
end