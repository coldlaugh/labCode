function [stress_force_yy]=stress(X,X0,B,S2,S,W)
    stress_force_yy = 0;
    BC = @(x) [x(1)-round(x(1)/W)*W,x(2)];
    for i = 1:size(B,1)
        n1 = B(i,1); n2 = B(i,2);
        if (any(S2==n1) && any(S==n2)) || (any(S==n1) && any(S2==n2))
            dr = BC(X(n1,:) - X(n2,:));
            r1 = norm(dr);
            r0 = norm(BC(X0(n1,:) - X0(n2,:)));
            stress_force_yy = stress_force_yy + (r1 - r0) / r1 * dr(2) / W;
        end
    end
end