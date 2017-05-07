function TopologicalLattice(domain_wall, notch_length, wall_distance)

%{

Lattice fracturing simulation.
Quasi static uniaxial pulling.

domain_wall = true / false
notch_length is the length of the notch
wall_distance is the domain wall distance
label is a unique label for file output


%}
%% declear global variables

%% coder extrinsic

% declear external matlab function for coder
coder.extrinsic('sprintf')
coder.extrinsic('initLattice')
% coder.extrinsic('breaking_threshold')
coder.extrinsic('break_bonds')
coder.extrinsic('figure_putput')
coder.extrinsic('stress')
coder.extrinsic('notch')
coder.extrinsic('mkdir')

% pre assign matrices with variable shapes %%% Now I need a placeholder
% function that create the correct dim for global variables.

% X =zeros(72000,2);  coder.varsize('X', [72000, 2], [1 0]);
% X0 = zeros(7200,2); coder.varsize('X0', [72000, 2], [1 0]);
% B = zeros(72000,2); coder.varsize('B', [72000, 2], [1 0]);
% S1 = zeros(7200,1); coder.varsize('S1', [7200, 1], [1 0]);
% S2 = zeros(7200,1); coder.varsize('S2', [7200, 1], [1 0]);
% S = zeros(1,72000); coder.varsize('S', [1, 72000], [0 1]);
% B_broken = B; coder.varsize('B_broken', [72000, 2], [1 0]);
% B_broken0 = B; coder.varsize('B_broken0', [72000, 2], [1 0]);
% G = zeros(72000,1); coder.varsize('G', [72000, 1], [1 0]);

% pre assign useful variables
W = 0.0;    %lattice width
H = 0.0;    %lattice height
snapshot = 0;   %snapshot counting
stress_force_yy = 0.0;  %yy component of stress matrix
avalanche = 0; %% avalanche counting
fileID = 0; 
t=0;



%% Management of DIR

% cd '/Users/leyou/Box Sync/MATLAB/TL/Fracture'
% set file name 
rng('shuffle');
label = randi(100000000);
if domain_wall
    mkdir('result/notch');
	file_avalanche = sprintf('result/notch/avalanche_wall_notch%d_walldistance_%d_label%d.txt',notch_length,wall_distance,label);
	file_stress = sprintf('result/notch/stress_wall_notch%d_walldistance_%d_label%d.txt',notch_length,wall_distance,label);
else
    mkdir('result/notch');
	file_avalanche = sprintf('result/notch/avalanche_control_notch%d_label%d.txt',notch_length,label);
	file_stress = sprintf('result/notch/stress_control_notch%d_label%d.txt',notch_length,label);
end

%% Main flow

% 1. initialize the lattice. 

show_plot = false; % whether to plot lattice while simulation is going on
B_broken = zeros(0,2); % record broken bonds
B_broken0 = zeros(0,2); % backup B_broken
% domain_wall = true;

latticeLength = 30;
latticeNumStack = 30;
%first give place holder to X,B,S1,S2,S,W,H
[X,B,S1,S2,S,lRest] = initPlaceHolder(latticeLength, latticeNumStack);

[X,B,S1,S2,S,W,H] = initLattice(latticeLength,latticeNumStack,show_plot,domain_wall,wall_distance);

[lRest] = bondLength(X,B,W);

X0 = X;  %% this line of code may be useless


% 1.5 Add notch (vertical, in the center)

%     [B] = notch(X,B,W,H,notch_length); % let's first ignore the notch
    
% 2. Assign random breaking threshold

    [G] = breaking_threshold(0.05,0.005,B);

% 3. Loop of strain

for strain = [0.00:0.001:0.2]
    a = sprintf('strain=%f',strain);
    a
%     strain
    snapshot = snapshot + 1;

% 4. Energy minimization

    flag_relaxation = true;
    B_broken0 = B_broken;
    while flag_relaxation 
		for i = [1:numel(S2)]
			j = S2(i);
        	X(j,2) = X0(j,2) + H * strain;
        end
        
        [X] = relaxation(X,B,S1,S2,lRest,W); % relaxation of position X, global variable in this function


% 5. Breaking of bonds
    
        [B,B_broken,G,flag_relaxation] = break_bonds(X,X0,B,B_broken,S,G,W);
    end
    

% 6. (Optional) Figure output
    
    if show_plot
        figure_output()
    end

% 7. Store results in files
    [stress_force_yy] = stress(X,X0,B,S2,S,W);
    avalanche = size(B_broken,1) - size(B_broken0,1);
    fileID = fopen(file_stress,'a');
    fprintf(fileID,'%16.8f %16.8f\n',strain,stress_force_yy);
    fclose(fileID);
    if avalanche > 0
        fileID = fopen(file_avalanche,'a');
        fprintf(fileID,'%16.8f %16.8f\n',strain,avalanche);
        fclose(fileID);
    end
    
end


end
%% Functions
