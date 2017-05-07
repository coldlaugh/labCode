function [X,B,S1,S2,S,R0]=initPlaceHolder(L,nstack) % place hold to creat the right dim for global variables
    % first layer:
    
        layer1=zeros(2*L,2);
    % second layer:
        layer2 = zeros(L,2);

    % third layer
        layer3=zeros(2*L,2);

    % forth layer
        layer4 = zeros(L,2);
        
        % list of bonds, site positions and C matrix 

    X=zeros(6*L*nstack,2);
    B = zeros(12*L*nstack-2*L,2);
    R0 = zeros(12*L*nstack-2*L,1);

    % site list
   

    % top and bottom boundary

    n_top = [size(X,1)-size(layer4,1)-size(layer3,1)+1:size(X,1)];
    n_bottom = [1:size(layer1,1)];

    % inside sites

    n_in = [size(layer1,1)+1:size(X,1)-size(layer4,1)-size(layer3,1)];
    
    S1 = n_bottom;
    S2 = n_top;
% 
%     S2 = sort(unique(S2));
%     S = setdiff([1:size(X,1)],S1);
%     S = setdiff(S,S2);
    S = n_in;
    
    
end