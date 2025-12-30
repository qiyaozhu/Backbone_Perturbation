% Calculate Lennard-Jones repulsive
% Glycine each residue has 7 atoms: N, H, CA, 1HA, 2HA, C, O
function fa_rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well)

fa_rep = 0;
N = size(coordinates,1); % number of atoms
n = N/7; % number of residues
threshold = 3.4473*2*sqrt(3);

% Pick Ca as the residue centers and calcuate their Manhattan distances 
for r1 = 1 : n-1
    for r2 = r1+1 : n
        % Compute atom-pair interaction for residues close enough
        if sum(abs(coordinates(7*r1-4,:) - coordinates(7*r2-4,:))) <= threshold
            for i = 1 : 7
                for j = 1 : 7
                    a1 = (r1-1)*7 + i;
                    a2 = (r2-1)*7 + j;

                    k = D(a1, a2);
                    
                    if k >= 4
                        % set the connectivity weight
                        if k == 4
                            w = 0.2;
                        else
                            w = 1;
                        end
                        
                        E_fa_rep = atom_pair_rep(norm(coordinates(a1,:) - coordinates(a2,:)), ...
                            LJ_radius(a1), LJ_well(a1), LJ_radius(a2), LJ_well(a2));
                        fa_rep = fa_rep + w*E_fa_rep;
                    end
                end
            end

        end
    end
end

% In the same residue, only H and O are separated by 4 bonds
w = 0.2;
for i = 1 : 7
    a1 = i*7-5;
    a2 = i*7;
    E_fa_rep = atom_pair_rep(norm(coordinates(a1,:) - coordinates(a2,:)), ...
                            LJ_radius(a1), LJ_well(a1), LJ_radius(a2), LJ_well(a2));
    fa_rep = fa_rep + w*E_fa_rep;
end
end


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