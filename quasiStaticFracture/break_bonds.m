function [B,B_broken,relaxation]=break_bonds(X,lRest,B,B_broken,S,G,W)
    
    coder.extrinsic('datasample');
    breakIndex = 0;
    
    lCurrent = bondLength(X,B,W);
    r = zeros(size(B,1),1);
    for i = 1:size(B,1)
        n1 = B(i,1); n2 = B(i,2);
        if n1>=0 && n2 >=0
            if any(S==n1) || any(S==n2)
                r(i) = lCurrent(i) - lRest(i);
            end
        end
    end

%     [maxforce maxforcebond] = max(r./G);
    forceRatio = r ./ G ./ lRest; % the ratio between the relative extension and the threshold
    relaxation = any(forceRatio > 1);
    maxRatio = max(forceRatio)
    if relaxation
        breakIndex = datasample(find(forceRatio),1); % randomly sample one bond that is about to break
        B_broken = [B_broken;
            B(breakIndex,:)];
        B(breakIndex,[1,2]) = [-1,-1]; % mark the bond as broken
    end  
end