%% input : bond list: B, bond position X, initial position X0, top boundary list S1, lower boundary list S2,strain 
% output: updated bond position X1,

function [X] = relaxation(X,B,S1,S2,lRest,W)


ns = size(X,1);
nb = size(B,1);
% S = setdiff([1:ns],S1,'sorted');
% S = setdiff(S,S2,'sorted'); % inside sites

Grad = zeros(size(X));

% rest length %%% change it to a input value, so we don't have to calculate
% it every time.

% r0 = zeros(nb,1);
% 
% for i = 1:nb
%     n1 = B(i,1); n2 = B(i,2);
%     dr = X0(n1,:) - X0(n2,:);
%     dr(1) = dr(1) - round(dr(1)/W)*W;
%     r0(i) = norm(dr);
% end

% velocity

vel = zeros(size(X));

% coefficients
coeff = 0.02;
f = coeff;
nstep = 0;
dt = 1e-4;
eps = 1e-8;  % the cut off of grad magnitude

while true
    
    % Optimization 1
    
    
    X = X + dt * vel - 0.5 * Grad * dt.^2;
    vel = vel - 0.5 * dt * Grad;
    
    
    % Gradient
    
    Grad = zeros(size(Grad));
    for i = 1:nb
        n1 = B(i,1); n2 = B(i,2);
        if n1 >= 0 && n2 >= 0
            dr = X(n1,:) - X(n2,:); 
            dr(1) = dr(1) - round(dr(1)/W)*W; % the vector of the bond
            r = norm(dr);  % the current bond length
            g = (r - lRest(i))/r * dr / lRest(i);  %spring constant is inverse perportional to rest length
            Grad(n1,:) = Grad(n1,:) + g;
            Grad(n2,:) = Grad(n2,:) - g;
        end
    end

    Grad(S1,:) = 0; Grad(S2,:) = 0;
    vel(S1,:) = 0; vel(S2,:) = 0;
    % Optimization 2

    vel = vel - 0.5 * dt * Grad;
    
    % Conditions
    
%     if mod(i,100) == 0
%         maxGrad = max(max(abs(Grad)))
%         maxvel = max(max(abs(vel)))
%     end
    
    maxGrad = max(max(abs(Grad)));
    
    if maxGrad <= eps
        break;
    end
    
    vf = - sum(sum(Grad .* vel));
    vv = sum(sum(vel .* vel));
    ff = sum(sum(Grad .* Grad));
    
    if vf < 0
        vel = zeros(size(vel));
        dt = dt * 0.5;
        f = coeff;
        nstep = 0;
    else
        nstep = nstep + 1;
        if nstep > 5
            dt = min([dt*1.1,5e-3]);
            f = coeff * 0.99;
        end
    end
    
    vel = vel*(1-f) - f * (vv/ff).^0.5 * Grad;
    
end


end