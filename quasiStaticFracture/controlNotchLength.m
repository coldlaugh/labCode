function mass_produce_notch(domain_wall, notch_length, wall_distance, label) %#codegen

%% coder extrinsic

coder.extrinsic('sprintf')
coder.extrinsic('initial')
coder.extrinsic('breaking_threshold')
coder.extrinsic('break_bonds')
coder.extrinsic('figure_putput')
coder.extrinsic('stress')
coder.extrinsic('notch')
X =zeros(72000,2);
coder.varsize('X', [72000, 2], [1 0]);
X0 = zeros(7200,2);
coder.varsize('X0', [72000, 2], [1 0]);
B = zeros(72000,2);
coder.varsize('B', [72000, 2], [1 0]);
S1 = zeros(7200,1);
coder.varsize('S1', [7200, 1], [1 0]);
S2 = zeros(7200,1);
coder.varsize('S2', [7200, 1], [1 0]);
S = zeros(1,72000);
coder.varsize('S', [1, 72000], [0 1]);
W = 0.0;
H = 0.0;
B_broken = B;
coder.varsize('B_broken', [72000, 2], [1 0]);
B_broken0 = B;
coder.varsize('B_broken0', [72000, 2], [1 0]);
G = zeros(72000,1);
coder.varsize('G', [72000, 1], [1 0]);
snapshot = 0;
stress_force_yy = 0.0;
avalanche = 0;
fileID = 0;
t=0;

%% massive production of fracture.

%% The first step is to sample around 20 configs for the two cases.
%  Thresholds will first be 0.05 +/- 0.01 for this moment.

%% Management of DIR

% cd '/Users/leyou/Box Sync/MATLAB/TL/Fracture'

if domain_wall
	file_avalanche = sprintf('data_notch/avalanche_wall_notch%d_walldistance_%dlabel%d.txt',notch_length,wall_distance,label);
	file_stress = sprintf('data_notch/stress_wall_notch%d_walldistance_%d_label%d.txt',notch_length,wall_distance,label);
else
	file_avalanche = sprintf('data_notch/avalanche_control_notch%d_label%d.txt',notch_length,label);
	file_stress = sprintf('data_notch/stress_control_notch%d_label%d.txt',notch_length,label);
end

%% Main flow

% 1. initialize the lattice. 

show_plot = false;
B_broken = zeros(0,2);
B_broken0 = zeros(0,2);
% domain_wall = true;
[X,B,S1,S2,S,W,H]=initial(50,80,show_plot,domain_wall,wall_distance);

X0 = X;


% 1.5 Add notch (vertical, in the center)

    [B] = notch(X,B,W,H,notch_length);
    
% 2. Assign random breaking threshold

    [G] = breaking_threshold(0.05,0.005,B,label);

% 3. Loop of strain

for strain = [0.00:0.001:0.2]
    
%     strain
    snapshot = snapshot + 1;

% 4. Energy minimization

    flag_relaxation = true;
    B_broken0 = B_broken;
    while flag_relaxation 
		for i = [1:numel(S2)]
			j = S2(i);
        	X(j,2) = X0(j,2) + X0(j,2) * strain;
		end
        [X] = relaxation(X,X0,B,S1,S2,W);


% 5. Breaking of bonds
    
        [B,B_broken,G,flag_relaxation]=break_bonds(X,X0,B,B_broken,S,G,W);
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
