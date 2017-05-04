%% input : bond list: B, bond position X, initial position X0, top boundary list S1, lower boundary list S2,strain 
% output: updated bond position X1,

function [X1] = relaxation(X,X0,B,S1,S2,W)

ns = size(X,1);
nb = size(B,1);
S = setdiff([1:ns],S1,'sorted');
S = setdiff(S,S2,'sorted'); % inside sites

Grad = zeros(size(X));

% rest length

r0 = zeros(nb,1);

for i = 1:nb
    n1 = B(i,1); n2 = B(i,2);
    dr = X0(n1,:) - X0(n2,:);
    dr(1) = dr(1) - round(dr(1)/W)*W;
    r0(i) = norm(dr);
end

% velocity

vel = zeros(size(X));

% coefficients
coeff = 0.5;
f = coeff;
nstep = 0;
dt = 1e-4;
eps = 1e-5;

while true
    
    % Optimization 1
    
    
    X = X + dt * vel - 0.5 * Grad * dt.^2;
    vel = vel - 0.5 * dt * Grad;
    
    
    % Gradient
    
    Grad = zeros(size(Grad));
    
    for i = 1:nb
        n1 = B(i,1); n2 = B(i,2);
        dr = X(n1,:) - X(n2,:);
        dr(1) = dr(1) - round(dr(1)/W)*W;
        r = norm(dr);
        Grad(n1,:) = Grad(n1,:) + (r - r0(i))/r * dr;
        Grad(n2,:) = Grad(n2,:) - (r - r0(i))/r * dr;
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
            nstep = nstep + 1;
        end
    end
    
    vel = vel*(1-f) - f * (vv/ff).^0.5 * Grad;
    
end

X1 = X;

end