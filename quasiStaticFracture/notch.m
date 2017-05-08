function [B] = notch(X,B,W,H,a) %create a notch in the lattice
    notch_center = [W/2 H/2];
    i = 1;
    for i = [1:size(B,1)]
        n1 = B(i,1); n2 = B(i,2);
        if n1 >= 0 && n2 >= 0   %
            if notch_boundary(X(n1,:),notch_center,a) && notch_boundary(X(n2,:),notch_center,a)
                B(i,:) = [-1,-1];
            end
        end
    end
end

function [e] = notch_boundary(x,c,a)
    xc = x - c;
    if abs(xc(1))<=2 && abs(xc(2))<(a/2)
        e = true;
    else
        e = false;
    end
    return;
end
    