clc;
clear all;

%%%%%%%%%%%%%%%%%% PARAMETERS AND PATHS TO BE CHANGED %%%%%%%%%%%%%%%%%%
addpath("top100fullStructure");
addpath("my_functions");

% pdb filename of the input structure
complex_name = "alpha_01_1a85_cyclic_001";

% Number of Monte Carlo trajectories for backbone perturbation
ip_size = 100;

% Output folder for the perturbed backbones
output_folder = "perturbed_backbones/";

% Target residues for protein binding site
target_res = [66,67,73,81,104,105,122,134];

% Chain ID for protein and peptide
peptide_chain = "A";
protein_chain = "B";

% Ideal CA-CA distance between the peptide and the protein target residues,
% rewards follows Gaussian distribution, mu and sigma are the mean and std
target_mu = 6;
target_sigma = 2;

% Run the perturbation function
perturbation(complex_name, target_res, ip_size, output_folder, protein_chain, peptide_chain, target_mu, target_sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function perturbation(complex_name, target_res, ip_size, output_folder, protein_chain, peptide_chain, target_mu, target_sigma)

% Thresholds for bond geometry deviation
length_thresh = 0.05; % Å
angle_thresh = 5; % degrees
omega_thresh = 10; % degrees

complex = pdb2mat(complex_name+".pdb");
X = complex.X;
Y = complex.Y;
Z = complex.Z;
atom_names = complex.atomName;
res_names = complex.resName;
res_num = complex.resNum;
element = complex.element;
chain_id = complex.chainID;

% Extract the peptide backbone
peptide_X = X(chain_id==peptide_chain);
peptide_Y = Y(chain_id==peptide_chain);
peptide_Z = Z(chain_id==peptide_chain);
peptide_coordinates = [peptide_X.', peptide_Y.', peptide_Z.'];
peptide_res_num = res_num(chain_id==peptide_chain);
peptide_atom_names = atom_names(chain_id==peptide_chain);

n = peptide_res_num(end);
peptide_backbone = zeros(n*4, 3);
for i = 1 : n
    peptide_backbone(i*4-3, :) = peptide_coordinates(peptide_res_num==i & peptide_atom_names=="N", :);
    peptide_backbone(i*4-2, :) = peptide_coordinates(peptide_res_num==i & peptide_atom_names=="CA", :);
    peptide_backbone(i*4-1, :) = peptide_coordinates(peptide_res_num==i & peptide_atom_names=="C", :);
    peptide_backbone(i*4, :) = peptide_coordinates(peptide_res_num==i & peptide_atom_names=="O", :);
end

% compute the phi, psi torsion angles
peptide_torsions = zeros(n*2, 1);
for i = 1 : n
    if i == 1
        C_pre = peptide_backbone(n*4-1, :);
    else
        C_pre = peptide_backbone(i*4-5, :);
    end
    N = peptide_backbone(i*4-3, :);
    CA = peptide_backbone(i*4-2, :);
    C = peptide_backbone(i*4-1, :);
    if i == n
        N_next = peptide_backbone(1, :);
    else
        N_next = peptide_backbone(i*4+1, :);
    end

    peptide_torsions(2*i-1) = torsion(C_pre, N, CA, C);
    peptide_torsions(2*i) = torsion(N, CA, C, N_next);
end

% Check bond geometry first
bond_check = true;

% compute the omega torsions, bond angles, and bond lengths
bond_lengths = zeros(n, 3); % N-CA, CA-C, C-N
bond_angles = zeros(n, 3); % C-N-CA, N-CA-C, CA-C-N
omega = zeros(n, 1);

% Compute omega angles
for i = 1 : n
    i_pre = i-1;
    i_pre = (i_pre==0)*n + (i_pre~=0)*i_pre;
    i_next = mod(i, n) + 1;

    C_pre = peptide_backbone(4*i_pre-1,:);
    N_curr  = peptide_backbone(4*i-3,:);
    CA_curr = peptide_backbone(4*i-2,:);
    C_curr  = peptide_backbone(4*i-1,:);
    N_next  = peptide_backbone(4*i_next-3,:);
    CA_next  = peptide_backbone(4*i_next-2,:);

    % Omega angle computation
    omega(i) = dihedral(CA_curr, C_curr, N_next, CA_next);

    % Bond lengths
    bond_lengths(i, 1) = norm(N_curr - CA_curr); % N-CA
    bond_lengths(i, 2) = norm(CA_curr - C_curr); % CA-C
    bond_lengths(i, 3) = norm(C_curr - N_next);  % C-N

    % Bond angles
    bond_angles(i, 1) = angle3(C_pre, N_curr, CA_curr);  % C-N-CA
    bond_angles(i, 2) = angle3(N_curr, CA_curr, C_curr);  % N-CA-C
    bond_angles(i, 3) = angle3(CA_curr, C_curr, N_next);  % CA-C-N
end

% Check for large omega deviations
omega_diff = abs(omega - 180);
omega_diff = min(omega_diff, 360-omega_diff);

for res = 1 : n
    if omega_diff(res) > omega_thresh
        bond_check = false;
        fprintf('Omega angle deviation: residue %d deviation %.3f degrees\n', res, omega_diff(res));
    end
end

% Canonical ideal bond lengths (Å) and angles (degrees)
ideal_lengths = [1.458, 1.524, 1.329]; % [N-CA, CA-C, C-N]
bond_type = ["N-CA", "CA-C", "C-N"];
ideal_angles  = [121.7, 111.2, 116.2]; % [C-N-CA, N-CA-C, CA-C-N]
angle_type = ["C-N-CA", "N-CA-C", "CA-C-N"];

% Compute deviations from ideal with non-canonical handling
length_deviation = bond_lengths-ideal_lengths;
angle_deviation  = bond_angles-ideal_angles;

for res = 1 : n
    for type = 1 : 3
        if length_deviation(res, type) > length_thresh
            bond_check = false;
            fprintf('Bond length deviation: residue %d bond %s deviation %.3f Å\n', res, bond_type(type), length_deviation(res, type));
        end
    end
end

for res = 1 : n
    for type = 1 : 3
        if angle_deviation(res, type) > angle_thresh
            bond_check = false;
            fprintf('Bond angle deviation: residue %d angle %s deviation %.3f degree\n', res, angle_type(type), angle_deviation(res, type));
        end
    end
end

% For peptides with good bond geometry, proceed to backbone perturbation
if bond_check
    % Build the scoring grid for rough computation of repulsive energy
    % between the peptide and the protein
    grid_dx = 1;
    threshold = ceil(4/grid_dx); % distance cutoff for repulsive energy
    padding = 12;

    minX = floor(min(X))-grid_dx*padding;
    maxX = ceil(max(X))+grid_dx*padding;
    nX = ceil((maxX-minX)/grid_dx);

    minY = floor(min(Y))-grid_dx*padding;
    maxY = ceil(max(Y))+grid_dx*padding;
    nY = ceil((maxY-minY)/grid_dx);

    minZ = floor(min(Z))-grid_dx*padding;
    maxZ = ceil(max(Z))+grid_dx*padding;
    nZ = ceil((maxZ-minZ)/grid_dx);

    % Extract protein coordinates
    protein_X = X(chain_id==protein_chain);
    protein_Y = Y(chain_id==protein_chain);
    protein_Z = Z(chain_id==protein_chain);
    protein_atom_names = atom_names(chain_id==protein_chain);
    protein_coordinates = [protein_X.', protein_Y.', protein_Z.'];
    protein_res_names = res_names(chain_id==protein_chain);
    protein_res_num = res_num(chain_id==protein_chain);
    protein_element = element(chain_id==protein_chain);

    % The neighborhood region that a grid point needs to check if protein atoms
    % reside in for computing repulsive energy
    neighbor = [];
    for i = -threshold : 1 : threshold
        for j = -threshold : 1 : threshold
            for k = -threshold : 1 : threshold
                pos = [i, j, k];
                if norm(pos) <= threshold
                    neighbor = [neighbor; pos];
                end
            end
        end
    end

    % LJ radius
    s_C = 2.0;
    s_N = 1.8;
    s_O = 1.5;
    s_S = 2.0;
    s_H = 1.0;

    % LJ well depth
    e_C = 0.06;
    e_N = 0.16;
    e_O = 0.16;
    e_S = 0.46;
    e_H = 0.02;

    % Repulsive energy weight for protein backbone and side-chain atoms
    w_bb = 1;
    w_sc = 0.2;

    % Record which grid cells contain protein atoms
    progrid = cell(nX, nY, nZ);
    for i = 1 : length(protein_X)
        xpos = floor((protein_X(i)-minX)/grid_dx)+1;
        ypos = floor((protein_Y(i)-minY)/grid_dx)+1;
        zpos = floor((protein_Z(i)-minZ)/grid_dx)+1;

        element = protein_element{i};
        if element == "C"
            atom.s = s_C;
            atom.e = e_C;
        elseif ismember(element, ["N", "N1+"])
            atom.s = s_N;
            atom.e = e_N;
        elseif ismember(element, ["O", "O1-"])
            atom.s = s_O;
            atom.e = e_O;
        elseif element == "S"
            atom.s = s_S;
            atom.e = e_S;
        elseif element == "H"
            atom.s = s_H;
            atom.e = e_H;
        else
            fprintf("No element found for "+element+"!\n");
        end

        name = protein_atom_names{i};
        if ismember(name, ["N", "H", "CA", "C", "O", "HA", "CB"])
            atom.w = w_bb;
        else
            atom.w = w_sc;
        end

        atom.pos = [protein_X(i), protein_Y(i), protein_Z(i)];
        progrid{xpos,ypos,zpos} = [progrid{xpos,ypos,zpos}, atom];
    end

    % Compute the repulsive energy for each grid cell, between the cell center
    % and all protein atoms in the neighborhood
    repulsive_grid_C = zeros(nX, nY, nZ);
    repulsive_grid_N = zeros(nX, nY, nZ);
    repulsive_grid_O = zeros(nX, nY, nZ);
    repulsive_grid_H = zeros(nX, nY, nZ);

    for x = threshold+1 : nX-threshold
        for y = threshold+1 : nY-threshold
            for z = threshold+1 : nZ-threshold
                % cell center coordinates
                center = [minX+(x-0.5)*grid_dx, minY+(y-0.5)*grid_dx, minZ+(z-0.5)*grid_dx];

                % for each neighbor cell, check the protein atoms reside in
                for nb = 1 : size(neighbor,1)
                    neighbor_pos = [x,y,z] + neighbor(nb, :);
                    atoms = progrid{neighbor_pos(1), neighbor_pos(2), neighbor_pos(3)};
                    for a = 1 : length(atoms)
                        atom = atoms(a);
                        repulsive_grid_C(x,y,z) = repulsive_grid_C(x,y,z) + atom.w * atom_pair_rep(norm(center-atom.pos), s_C, e_C, atom.s, atom.e);
                        repulsive_grid_N(x,y,z) = repulsive_grid_N(x,y,z) + atom.w * atom_pair_rep(norm(center-atom.pos), s_N, e_N, atom.s, atom.e);
                        repulsive_grid_O(x,y,z) = repulsive_grid_O(x,y,z) + atom.w * atom_pair_rep(norm(center-atom.pos), s_O, e_O, atom.s, atom.e);
                        repulsive_grid_H(x,y,z) = repulsive_grid_H(x,y,z) + atom.w * atom_pair_rep(norm(center-atom.pos), s_H, e_H, atom.s, atom.e);
                    end
                end
            end
        end
    end

    % CA atom coordinates of the target residues for specificity
    target = [];
    for i = 1 : length(protein_X)
        if ismember(protein_res_num(i),target_res) && protein_atom_names(i)=="CA"
            target = [target; protein_X(i), protein_Y(i), protein_Z(i)];
        end
    end

    % Negative scores for target region promotion
    % Use gaussian distribution for calculating score
    rewardgrid = zeros(nX, nY, nZ, size(target,1));
    target_grid = zeros(size(target));
    min_s = -10;

    for i = 1 : size(target,1)
        % Find the grid cells for the target residues
        xpos = floor((target(i,1)-minX)/grid_dx)+1;
        ypos = floor((target(i,2)-minY)/grid_dx)+1;
        zpos = floor((target(i,3)-minZ)/grid_dx)+1;
        target_grid(i,:) = [xpos, ypos, zpos];

        % Assign rewards to neighbors
        max_dist = target_mu + 3*target_sigma;
        [X,Y,Z]=ndgrid(-ceil(max_dist/grid_dx):ceil(max_dist/grid_dx), -ceil(max_dist/grid_dx):ceil(max_dist/grid_dx), -ceil(max_dist/grid_dx):ceil(max_dist/grid_dx));
        X = X(:);
        Y = Y(:);
        Z = Z(:);

        for k = 1 : length(X)
            if xpos+X(k)>=1 && xpos+X(k)<=nX && ypos+Y(k)>=1 && ypos+Y(k)<=nY && zpos+Z(k)>=1 && zpos+Z(k)<=nZ
                distance = norm([X(k), Y(k), Z(k)]) * grid_dx;
                rewardgrid(xpos+X(k), ypos+Y(k), zpos+Z(k), i) = normpdf(distance,target_mu,target_sigma) / normpdf(target_mu,target_mu,target_sigma) * min_s;
            end
        end
    end

    % Ramachandran plot for symmetric glycine
    Boundaries = load('ramabin_glycine.mat').Boundaries.';

    % atom properties
    [LJ_radius, LJ_well, ~, ~, ~, ~] = FA_parameter(n); % parameters
    D = n_bonds(n); % number of bonds between any atom pair

    % construct the peptide from the torsion angles
    C_start = peptide_backbone(end-1,:);
    N_start = peptide_backbone(1,:);
    CA_start = peptide_backbone(2,:);

    x = N_start - C_start;
    x = x/norm(x);
    u = CA_start - C_start;
    y = u - dot(u,x)*x;
    y = y/norm(y);
    z = cross(x, y);

    % % Sanity check for reconstructing the original input peptide
    % [coordinates, cyc] = peptide_cyc(peptide_torsions, C_start, N_start, x, y, z, bond_lengths, deg2rad(bond_angles), deg2rad(omega));
    % plot_complex(coordinates, protein_coordinates, protein_atom_names, protein_res_names, protein_res_num, protein_element, n, "test.pdb");

    % Simulated annealing parameters
    w_rep = 0.3;
    w_target = 0.3;
    w_hbd_ramping = 10;
    w_cyc = 12;

    t0_ramping = 50;
    k0 = 0.1;
    b = 18;
    c = 40;

    M = 10000;
    ramp0 = 1;

    % criteria for good candidates
    rep_cri = 3*n;
    target_cri = -7*size(target_grid,1);
    cyc_cri = 0.3;
    count_cri = ceil(n/3);

    % Compute the energies for the original input peptide
    [coordinates, cyc] = peptide_cyc(peptide_torsions, C_start, N_start, x, y, z, bond_lengths, deg2rad(bond_angles), deg2rad(omega));
    rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
    E_target_rep = target_rep(coordinates, repulsive_grid_C, repulsive_grid_N, repulsive_grid_O, repulsive_grid_H, minX, minY, minZ, grid_dx, nX, nY, nZ);
    rep_total = rep + E_target_rep;
    E_target_reward = target_reward(coordinates, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid);
    [hbond_ramping, count, oversat] = E_hbond_ramping(coordinates, ramp0);
    E_total_ramping = w_cyc*cyc + w_rep*rep_total + w_hbd_ramping*hbond_ramping;
    display("Original peptide total="+E_total_ramping+", cyc="+cyc+", peptide_rep="+rep+", target_rep="+E_target_rep+", target_reward="+E_target_reward+", hbond="+hbond_ramping+", hcount="+count);

    % Start Monte Carlo trajectories
    BEST_CAND_DATA = cell(1,ip_size);

    parfor repeat = 1 : ip_size
        best_cand = [];
        best_score = [];
        best_data = [];
        best_coordinates = [];

        E_total_ramping = 1000;

        accept_ramping = 0;
        angles = peptide_torsions;
        cand = 0;

        for i = 1 : M
            % generate new random move
            k = k0/(1+b*i/M);
            p = random_move(k,n);
            [p_new, angles_new] = feasible(p, angles, Boundaries);

            [coordinates, cyc] = peptide_cyc(angles_new, C_start, N_start, x, y, z, bond_lengths, deg2rad(bond_angles), deg2rad(omega));
            rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
            E_target_rep = target_rep(coordinates, repulsive_grid_C, repulsive_grid_N, repulsive_grid_O, repulsive_grid_H, minX, minY, minZ, grid_dx, nX, nY, nZ);
            rep_total = rep + E_target_rep;
            E_target_reward = target_reward(coordinates, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid);
            ramp = ramp0*(M-i)/M;
            [hbond_ramping, count, oversat] = E_hbond_ramping(coordinates, ramp);
            E_total_ramping_new = w_cyc*cyc + w_rep*rep_total + w_target*E_target_reward + w_hbd_ramping*hbond_ramping;

            pass = false;
            temp = t0_ramping/(1+c*i/M);
            if E_total_ramping_new <= E_total_ramping
                pass = true;
            else
                prob = exp(1)^((E_total_ramping-E_total_ramping_new)/temp);
                if rand <= prob
                    pass = true;
                end
            end

            % pass all layers, this random move is accepted
            if pass
                angles = angles_new;
                E_total_ramping = E_total_ramping_new;
                accept_ramping = accept_ramping + 1;

                if rep_total <= rep_cri && cyc <= cyc_cri && count >= count_cri && oversat <= 0 && E_target_reward <= target_cri
                    cand = cand + 1;
                    display("SA"+repeat+"_cand"+cand+", total="+E_total_ramping+", cyc="+cyc+", peptide_rep="+rep+", target_rep="+E_target_rep+", target_reward="+E_target_reward+", hbond="+hbond_ramping+", hcount="+count);
                    best_cand = [best_cand, cand];
                    best_score = [best_score, E_total_ramping];
                    best_data = [best_data; cyc, rep_total, count];
                    best_coordinates = [best_coordinates, coordinates];
                end
            end
        end

        if ~isempty(best_score)
            [~, best_ind] = mink(best_score, 1);
            BEST_CAND_DATA{repeat}.data = best_data(best_ind, :);
            BEST_CAND_DATA{repeat}.coordinates = best_coordinates(:, reshape([best_ind*3-2; best_ind*3-1; best_ind*3], [], 1));
        end

        display("SA = "+repeat+", Ramp accept = "+accept_ramping);
    end

    % Only take perturbed backbones that have 1~2 Å from the original
    total_cand = 0;
    coordinates = [];
    for repeat = 1 : ip_size
        if ~isempty(BEST_CAND_DATA{repeat})
            coor = BEST_CAND_DATA{repeat}.coordinates;
            for cand = 1 : size(coor,2)/3
                CO = coor(:,3*cand-2:3*cand);
                co = CO(repelem(1:7:7*n,4)+repmat([0,2,5,6],1,n),:);
                [~, ~, rmsd] = kabsch_algorithm(co, peptide_backbone);
                rmsd
                if rmsd <= 2
                    total_cand = total_cand + 1;
                    coordinates = [coordinates, CO];
                end
            end
        end
    end

    for cand_ind = 1 : total_cand
        filename = output_folder+complex_name+"_Perturb"+cand_ind+".pdb";
        plot_complex(coordinates(:, cand_ind*3-2:cand_ind*3), protein_coordinates, protein_atom_names, protein_res_names, protein_res_num, protein_element, n, filename, peptide_chain, protein_chain);
    end
end
end

% Helper function for computing the torsion angle between four points
function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
x = dot(n1, n2);
y = dot(cross(n1, n2), b2/norm(b2));
chi = atan2(y, x);
end


% Plot the peptide backbone
function plot_complex(peptide_coordinates, protein_coordinates, protein_atomName, protein_resName, protein_resNum, protein_element, n, filename, peptide_chain, protein_chain)
peptide.X = [peptide_coordinates(:,1).', protein_coordinates(:,1).'];
peptide.Y = [peptide_coordinates(:,2).', protein_coordinates(:,2).'];
peptide.Z = [peptide_coordinates(:,3).', protein_coordinates(:,3).'];

peptide.atomNum = 1 : length(peptide.X);
peptide.atomName = [repmat({"N", "H", "CA", "1HA", "2HA", "C", "O"}, 1, n), protein_atomName];

peptide.resName = [repmat({"GLY"}, 1, 7*n), protein_resName];
peptide.resNum = [repelem(1:n, 7), protein_resNum];
peptide.element = [repmat({"N", "H", "C", "H", "H", "C", "O"}, 1, n), protein_element];
peptide.chainID = [repmat({peptide_chain}, 1, size(peptide_coordinates,1)), repmat({protein_chain}, 1, length(protein_atomName))];

peptide.outfile = filename;
file = fopen(filename, "w");
fclose(file);
mat2pdb(peptide);
end


% Atom pair repulsive energy calculator
function E_fa_rep = atom_pair_rep(d, s_i, e_i, s_j, e_j)
epsilon = sqrt(e_i*e_j);
sigma = s_i + s_j;

if d <= 0.6*sigma
    m = 20*epsilon/sigma * (-(5/3)^12 + (5/3)^6);
    b = epsilon * (13*(5/3)^12 - 14*(5/3)^6 + 1);
    E_fa_rep = m*d+b;
elseif d <= sigma
    E_fa_rep = epsilon*((sigma/d)^12 - 2*(sigma/d)^6 + 1);
else
    E_fa_rep = 0;
end
end


% Function to calculate repulsive energy between peptide and protein
function score = target_rep(coordinates, repulsive_grid_C, repulsive_grid_N, repulsive_grid_O, repulsive_grid_H, minX, minY, minZ, grid_dx, nX, nY, nZ)
score = 0;

for i = 1 : size(coordinates,1)
    xpos = floor((coordinates(i,1)-minX)/grid_dx)+1;
    ypos = floor((coordinates(i,2)-minY)/grid_dx)+1;
    zpos = floor((coordinates(i,3)-minZ)/grid_dx)+1;

    % Depending on the atom type, different penalty grids, atom order in
    % backbone coordinates are N, H, CA, 1HA, 2HA, C, O
    if xpos>=1 && xpos<=nX && ypos>=1 && ypos<=nY && zpos>=1 && zpos<=nZ
        if mod(i, 7) == 3 || mod(i, 7) == 6
            score = score + repulsive_grid_C(xpos,ypos,zpos);
        elseif mod(i, 7) == 1
            score = score + repulsive_grid_N(xpos,ypos,zpos);
        elseif mod(i, 7) == 0
            score = score + repulsive_grid_O(xpos,ypos,zpos);
        else
            score = score + repulsive_grid_H(xpos,ypos,zpos);
        end
    end
end
end


% Function to sum rewards for all protein target sites.
% For each site, record the highest reward received from peptide CA atoms.
function score = target_reward(coordinates, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid)
n = size(coordinates,1) / 7;
target_score = zeros(size(target_grid,1),1);

% Only CA atoms contribute to rewards
for i = 1 : n
    xpos = floor((coordinates(i*7-4,1)-minX)/grid_dx)+1;
    ypos = floor((coordinates(i*7-4,2)-minY)/grid_dx)+1;
    zpos = floor((coordinates(i*7-4,3)-minZ)/grid_dx)+1;

    if xpos>=1 && xpos<=nX && ypos>=1 && ypos<=nY && zpos>=1 && zpos<=nZ
        rewards = rewardgrid(xpos,ypos,zpos,:);
        for t = 1 : length(target_score)
            if rewards(t) < target_score(t)
                target_score(t) = rewards(t);
            end
        end
    end
end

score = sum(target_score);
end


function [R, t, rmsd] = kabsch_algorithm(P, Q)

% Calculate the centroids of the two sets of points
centroid_P = mean(P, 1);
centroid_Q = mean(Q, 1);

% Center the points by subtracting their centroids
P_centered = P - centroid_P;
Q_centered = Q - centroid_Q;

% Calculate the covariance matrix of the centered points
covariance_matrix = P_centered' * Q_centered;

% Calculate the optimal rotation matrix using singular value decomposition (SVD)
[U, ~, V] = svd(covariance_matrix);
rotation_matrix = V * U';

% If the determinant of the rotation matrix is negative, we need to flip one axis
if det(rotation_matrix) < 0
    V(:, 3) = -V(:, 3);
    rotation_matrix = V * U';
end

% Calculate the translation vector
translation_vector = centroid_Q' - rotation_matrix * centroid_P';

% Apply the rotation and translation to the original set of points
P_aligned = (rotation_matrix * P')' + translation_vector';

% Calculate the root-mean-square deviation (RMSD) between the aligned points
rmsd = sqrt(sum(sum((Q - P_aligned).^2)) / size(P, 1));

% Output the rotation matrix, translation vector, and RMSD
R = rotation_matrix;
t = translation_vector;
end


% Helper function
function angle = dihedral(p1, p2, p3, p4)
b1 = p2 - p1;
b2 = p3 - p2;
b3 = p4 - p3;

n1 = cross(b1, b2);
n2 = cross(b2, b3);
m1 = cross(n1, b2/norm(b2));

x = dot(n1, n2);
y = dot(m1, n2);

angle = -atan2d(y, x);
end


function ang = angle3(a,b,c)
ba = a - b;
bc = c - b;
ang = acosd(dot(ba, bc) / (norm(ba)*norm(bc)));
end