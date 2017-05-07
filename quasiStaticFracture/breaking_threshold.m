function [G] = breaking_threshold(mu,sigma,B)
    rng('shuffle')
    G = normrnd(mu,sigma,size(B,1),1);
end