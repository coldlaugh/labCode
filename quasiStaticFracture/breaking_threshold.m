function [G] = breaking_threshold(mu,sigma,B,seed)
    t = datetime('now');
    e = 43 * minute(t) + 31 * second(t) + 21*seed;
    for i = 1:ceil(e)
        normrnd(mu,sigma,29,1);
    end
    G = normrnd(mu,sigma,size(B,1),1);
end