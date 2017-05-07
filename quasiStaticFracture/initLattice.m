function [X,B,S1,S2,S,W,H]=initLattice(L,nstack,show_plot,domain_wall,wall_distance)


    l = 2*[0.3 0.7 0.56; 0.3 0.7 0.56];  %% unit cell bond length
    theta = pi+0.2;  % angle between two triangles
    
    if domain_wall
        domain = @(x) (mod(x,wall_distance)<(wall_distance/2)); % define polarization direction
    else
        domain = @(x) true;
    end
    
    regular_kagome = false;  % if the lattice is a regular kagome
    if regular_kagome
        l = 2*[1 1 1; 1 1 1]; 
        theta = pi*2/3;
    end
    
    
    phi = acos((l(:,1).^2 + l(:,2).^2 - l(:,3).^2)./(2*l(:,1).*l(:,2)));  % the angle of the triangles.
    upper = [0 0;l(1,1) 0;cos(phi(1))*l(1,2) sin(phi(1))*l(1,2)];  % position of sites for upper triangle
    lower = [0 0;l(2,1) 0;cos(phi(2))*l(2,2) -sin(phi(2))*l(2,2)]; % lower triangle

    upper = upper * [cos(theta/2) sin(theta/2);-sin(theta/2) cos(theta/2)]; % rotate upper triangle by theta/2
    lower = lower * [cos(theta/2) -sin(theta/2);sin(theta/2) cos(theta/2)]; % rotate upper triangle by -theta/2

    upper2 = [-upper(:,1)+min(upper(:,1))*0.5 upper(:,2)];  % unit cell for another domain
    lower2 = [-lower(:,1)+min(lower(:,1))*0.5 lower(:,2)];
    
    % Lattice vectors

    T1 = lower(2,:) - upper(3,:); 
    T2 = upper(2,:) - lower(3,:);
    T3 = T1 + T2;
    
    
    % for a ref, see image resources. "unitCell.eps"
    
    % keep an index for layers
    n = 1;
    % first layer:
    
        layer1=[];
        layerIndex1 = [];
        for i = 1 : L
            if domain(i)
                layer1 = [layer1;lower([3 2],1)+T3(1)*i lower([3 2],2)+T3(2)*i];
            else
                layer1 = [layer1;lower2([2 3],1)+T3(1)*i lower2([2 3],2)+T3(2)*i];
            end
            layerIndex1 = [layerIndex1;n];
            n = n+1;
            layerIndex1 = [layerIndex1;n];
            n = n+1;
        end
    
    
    % second layer:
        layer2 = [];
        layerIndex2 = [];
        for i = 1 : L
            if domain(i)
                layer2 = [layer2;lower(1,1)+T3(1)*i lower(1,2)+T3(2)*i];
            else
                layer2 = [layer2;lower2(1,1)+T3(1)*i lower2(1,2)+T3(2)*i];
            end
            layerIndex2 = [layerIndex2;n];
            n = n + 1;
        end

    % third layer
        layer3=[];
        layerIndex3 = [];
        flag = true;
        for i = 1 : L
            if domain(i)
                layer3 = [layer3;
                    [upper(3,1) upper(3,2)]+T3*i;
                    [upper(2,1) upper(2,2)]+T3*i;
                    ];
            else
                layer3 = [layer3;
                    [upper2(2,1) upper2(2,2)]+T3*(i);
                    [upper2(3,1) upper2(3,2)]+T3*(i);
                    ];
            end
            layerIndex3 = [layerIndex3;n];
            n = n + 1;
            layerIndex3 = [layerIndex3;n];
            n = n + 1;
        end


    % forth layer
        layer4 = [];
        layerIndex4 = [];
        for i = 1 : L
            if domain(i)
                layer4 = [layer4;lower(1,1)+T3(1)*i-T1(1) lower(1,2)+T3(2)*i-T1(2)];
            else
                layer4 = [layer4;lower2(1,1)+T3(1)*(i-1)+T2(1) lower2(1,2)+T3(2)*(i-1)+T2(2)];  % this (i-1) could be ajusted.
            end
            layerIndex4 = [layerIndex4;n];
            n = n + 1;
        end
        
        % list of bonds, site positions and C matrix 

    B = [];
    B_boundary = [];
    B_h_boundary = [];
    X = [layer1;layer2;layer3;layer4];
    Ns = size(X,1);  % number of sites

    T = (-T1 + T2);
     % create bond list, the lattice has nstack*4 layers
    for k = 0 : nstack-1
        
        for i = 2:size(layer1,1)
%             B = [B;i-1+k*Ns i+k*Ns];
            B = [B; layerIndex1(i-1)+k*Ns, layerIndex1(i)+k*Ns];
        end


        % periodic BC

%         B = [B; i+k*Ns 1+k*Ns ];
        B = [B; layerIndex1(end)+k*Ns, layerIndex1(1)+k*Ns];



        for i = 1:size(layer1,1)
%             B = [B;i+k*Ns size(layer1,1)+ceil(i/2)+k*Ns];
            B = [B;layerIndex1(i)+k*Ns, layerIndex2(ceil(i/2))+k*Ns];
        end
        
        
        %PBC ?
        
        
        
        for i = 2:size(layer1,1)
%             B = [B;size(layer1,1)+size(layer2,1)+i-1+k*Ns size(layer1,1)+size(layer2,1)+i+k*Ns];

            B = [B; layerIndex3(i-1)+k*Ns, layerIndex3(i)+k*Ns];
        end


        % periodic BC

%         B = [B;size(layer1,1)+size(layer2,1)+i+k*Ns size(layer1,1)+size(layer2,1)+1+k*Ns];
        B = [B; layerIndex3(end)+k*Ns,layerIndex3(1)+k*Ns];


        for i = 1:size(layer3,1)
%             B = [B;size(layer1,1)+size(layer2,1)+i+k*Ns size(layer1,1)+ceil(i/2)+k*Ns];
            B = [B; layerIndex3(i)+k*Ns, layerIndex2(ceil(i/2))+k*Ns];
        end

        for i = 1:size(layer3,1)-1
%             B = [B;size(layer1,1)+size(layer2,1)+i+k*Ns size(layer1,1)+size(layer2,1)+size(layer3,1)+ceil((i+1)/2)+k*Ns];
            B = [B; layerIndex3(i)+k*Ns, layerIndex4(ceil((i+1)/2))+k*Ns];
        end
        
        % periodic BC
        
%         B = [B;size(layer1,1)+size(layer2,1)+size(layer3,1)+k*Ns size(layer1,1)+size(layer2,1)+size(layer3,1)+ceil((1+1)/2)+k*Ns];
        
        B = [B; layerIndex3(end)+k*Ns,layerIndex4(1)+k*Ns];

        if k ~= (nstack - 1)
            for i = 1:size(layer1,1)-1
%                 B = [B; i+Ns+k*Ns size(layer1,1)+size(layer2,1)+size(layer3,1)+ceil((i+1)/2)+k*Ns];
                B = [B; layerIndex1(i)+(k+1)*Ns, layerIndex4(ceil((i+1)/2))+k*Ns];
            end

            % periodic BC

%             B = [B; size(layer1,1)+Ns+k*Ns size(layer1,1)+size(layer2,1)+size(layer3,1)+1+k*Ns];
            B = [B; layerIndex1(end)+(k+1)*Ns, layerIndex4(1)+k*Ns];
            

        end
    end
    
    % site list
    X=[];
    for k = 0 : nstack-1
        X = [X;
            [layer1;layer2;layer3;layer4] + k*repmat(T,[Ns 1])];
    end

    % top and bottom boundary

    n_top = [size(X,1)-size(layer4,1)-size(layer3,1):size(X,1)];
    n_bottom = [1:size(layer1,1)];

    % inside sites

    n_in = [size(layer1,1)+1:size(X,1)-size(layer4,1)-size(layer3,1)-1];
    
%     S1 = [1:size([layer1;layer2;layer3],1)]';
    S1 = n_bottom;
%     S2 = [size(X,1)-size([layer1;layer2;layer3;layer4],1)+1:size(X,1)]';
    S2 = n_top;
%     H = mean(X(S2,2))-mean(X(S1,2));
%     W = (max(X(S1,1))-min(X(S1,1))) + (X(2,1)-X(1,1));
    H = T2(2) * (nstack+1); % the height of the box. No PBC in this direction.
    W = T3(1) * (L+1); % the width of the periodic box

    S2 = sort(unique(S2));
    S = setdiff([1:size(X,1)],S1);
    S = setdiff(S,S2);
    
    
    if show_plot
        % figure setting
        f=figure('Units','inches','Position',[0.2 0.2 L*5/10 8]);hold on; axis equal;
        haxis = gca;
        axis equal;
        BoxPos = [-1 1 L*4/10 6];
        set(haxis, ...
            'Units'              , 'inches',...
            'Position'           , BoxPos,...
            'XTick'              , [],...
            'YTick'              , [],...
            'XLim'               , [X(1,1),W-X(1,1)],...
            'YLim'               , [X(1,2),1.6*H],...
            'Box'                , 'on',...
            'LineWidth'          , 1);
        set(gcf,'color','white');

        curve = axes('Parent',f);
        set(curve, ...
            'Units'              , 'inches',...
            'Position'           , [BoxPos(3)-0.8 1 BoxPos(3)/10 BoxPos(4)],...
            'YLim'               , [0,0.5],...
            'Box'                , 'on',...,
            'LineWidth'          , 1);

%         if (cond == 1) || (cond == 2) || (cond == 3)
            set(curve, ...
                'XLim'               , [0 0.01]);
%         end

        xlabel(curve,'stress','FontSize',20,'FontName','Arial');
        ylabel(curve,'strain','FontSize',20,'FontName','Arial');
    end
    
end