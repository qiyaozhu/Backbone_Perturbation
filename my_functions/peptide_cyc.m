%% Glycine each residue has 7 atoms: N, H, CA, 1HA, 2HA, C, O
% @ angles: an input vector containing phi_1, psi_1, phi_2, psi_2, ..., phi_n, psi_n
% @ C_start, N_start: coordinates of the two starting atoms
% @ x, y, z: the coordinate system of the peptide
% @ bond lengths: n rows of [N-CA, CA-C, C-N] lengths
% @ bond angles: n rows of [C-N-CA, N-CA-C, CA-C-N] angles
% @ omega: the n omega torsion angles
% @ coordinates: the output atom coordinate matrix of size 7n x 3
% @ error: the cyclic error measuring how far this backbone is away from a perfect closure
% Note: To compute the cyclic error, we add a virtual atom C before the N
% terminus, and two virtual atoms N and CA after the C terminus.

function [coordinates, error] = peptide_cyc(angles, C_start, N_start, x, y, z, bond_lengths, bond_angles, omega)
% bond angles for adding H and O atoms
COR = deg2rad(120.8);
NHR = deg2rad(119.2);
NCaHa = deg2rad(109.5);

% bond lengths for adding H and O atoms
CO = 1.231;
NH = 1.010;
CaHa = 1.090;

% dihedral angles for 1HA and 2HA atoms on CA
CNCaH = deg2rad(121.4);
HNCaH = deg2rad(117.2);

% pre-computed translation vectors for adding the H atoms
NH_trans = [NH*cos(pi-NHR); -NH*sin(pi-NHR); 0];

% initialize variables
n = length(angles) / 2;
coordinates = zeros(7*n, 3);
coordinates(1,:) = N_start;
N = N_start;
rotation = eye(3);

% get the backbone coordinates
for i = 1 : n
    t = rotation*NH_trans;
    H = N + t(1)*x + t(2)*y + t(3)*z;
    
    rotation = rotation*T(bond_angles(i,1))*R(angles(2*i-1));
    t = rotation*[bond_lengths(i,1); 0; 0];
    Ca = N + t(1)*x + t(2)*y + t(3)*z;

    rotation_O = rotation;
    rotation = rotation*T(bond_angles(i,2))*R(angles(2*i));
    t = rotation*[bond_lengths(i,2); 0; 0];
    C = Ca + t(1)*x + t(2)*y + t(3)*z;

    t = rotation_O*T(bond_angles(i,2))*R(angles(2*i)-pi)*T(COR)*[CO; 0; 0];
    O = C + t(1)*x + t(2)*y + t(3)*z;

    rotation = rotation*T(bond_angles(i,3))*R(omega(i));
    t = rotation*[bond_lengths(i,3); 0; 0];
    N = C + t(1)*x + t(2)*y + t(3)*z;

    coordinates(7*i-5,:) = H;
    coordinates(7*i-4,:) = Ca;
    coordinates(7*i-1,:) = C;
    coordinates(7*i,:) = O;

    % get 1HA coordinates from C-N-CA-1HA
    coordinates(7*i-3,:) = kinematic_chain(C, coordinates(7*i-6,:), Ca, CaHa, NCaHa, CNCaH);

    % get 2HA coordinates from 1HA-N-CA-2HA
    coordinates(7*i-2,:) = kinematic_chain(coordinates(7*i-3,:), coordinates(7*i-6,:), Ca, CaHa, NCaHa, HNCaH);

    if i == n
        % Fix the last O position due to imperfect closure
        last_psi = torsion(coordinates(7*n-6,:), Ca, C, N_start);
        t = rotation_O*T(bond_angles(i,2))*R(last_psi-pi)*T(COR)*[CO; 0; 0];
        coordinates(7*n,:) = C + t(1)*x + t(2)*y + t(3)*z;

        % Compute Virtual_N and Virtual_CA
        Virtual_N = N;
        rotation = rotation*T(bond_angles(1,1));
        t = rotation*[bond_lengths(1,1); 0; 0];
        Virtual_CA = Virtual_N + t(1)*x + t(2)*y + t(3)*z;
    else
        coordinates(7*i+1,:) = N;
    end
end

% To fix the first H, need to go in reverse order, C1-CA1-N1-C_last
x = coordinates(3,:)-coordinates(6,:);
x = x/norm(x);
u = coordinates(1,:) - coordinates(6,:);
y = u - dot(u,x)*x;
y = y/norm(y);
z = cross(x, y);

phi_first = torsion(coordinates(6,:), coordinates(3,:), coordinates(1,:), coordinates(7*n-1,:));
t = T(bond_angles(1,2))*R(phi_first)*NH_trans;
H = coordinates(1,:) + t(1)*x + t(2)*y + t(3)*z;
coordinates(2,:) = H;

% Compute the cyclic error
f = @(x) (x<=1)*x + (x>1)*(x^2);
error = f(norm(coordinates(7*n-1,:)-C_start)) + f(norm(Virtual_N-N_start)) + f(norm(Virtual_CA-coordinates(3,:)));
end


% Helper function for computing rotation matrix caused by bond angles
function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end

% Helper function for computing rotation matrix caused by torsion angles
function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
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

% For a chain of four atoms p1-p2-p3-p4, compute the coordinates of p4,
% given: the coordinates of p1, p2, p3,
%        the bond length l of p3-p4,
%        the bond angle ba of p2-p3-p4,
%        the dihedral angle da of p1-p2-p3-p4.
function p4 = kinematic_chain(p1, p2, p3, l, ba, da)
x = p2 - p1;
x = x / norm(x);
u = p3 - p1;
y = u - dot(u,x)*x;
y = y / norm(y);
z = cross(x, y);

theta = acos(dot(p1-p2,p3-p2) / (norm(p1-p2)*norm(p3-p2)));

t = T(theta) * R(da) * T(ba) * [l; 0; 0];
p4 = p3 + t(1)*x + t(2)*y + t(3)*z;
end