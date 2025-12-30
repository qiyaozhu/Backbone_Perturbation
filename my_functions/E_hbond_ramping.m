%% This commented out section is for visualization of the energy functions used to compute hydrogen bond energies
% clc;
% clear all;
% 
% fs = 20;
% HA_max = 6;
% 
% xmin_HA = 1.384038127;
% xmax_HA = 2.998103943;
% fmin_HA = 1.1;
% fmax_HA = 1.1;
% C_HA = [-0.5307601, 6.47949946, -22.39522814, -55.14303544, 708.3094524, -2619.493182, 5227.88058, -6043.312116, 3806.046762, -1007.660241];
% xroot_HA = 2;
% 
% d_HA = linspace(0,HA_max,1000);
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% for r = 1 : length(ramp)
%     E_HA = zeros(1,length(d_HA));
%     for i = 1 : length(d_HA)
%         E_HA(i) = polynomial_HA(d_HA(i), xmin_HA, xroot_HA, xmax_HA, HA_max, fmin_HA, fmax_HA, C_HA, ramp(r));
%     end
%     plot(d_HA, E_HA, "LineWidth", 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Distance HA");
% ylabel("E_{HA}");
% fontsize(fs,"points");
% 
% 
% xmin_AHD = 1.143564639;
% xmax_AHD = 3.1416;
% fmin_AHD = 1.1;
% fmax_AHD = 1.1;
% C_AHD = [0.47683259, -9.54524724, 83.62557693, -420.5586777, 1337.193549, -2786.262657, 3803.178227, -3278.628799, 1619.041162, -347.5015791];
% 
% AHD = linspace(0,pi,1000);
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% for r = 1 : length(ramp)
%     E_AHD = zeros(1,length(d_HA));
%     for i = 1 : length(d_HA)
%         E_AHD(i) = polynomial_AHD(AHD(i), xmin_AHD, xmax_AHD, fmin_AHD, fmax_AHD, C_AHD, ramp(r));
%     end
%     plot(AHD, E_AHD, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Bond Angle AHD");
% ylabel("E_{AHD}");
% fontsize(fs,"points");
% 
% BAH = linspace(0,pi,1000);
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% for r = 1 : length(ramp)
%     E_BAH = zeros(1,length(BAH));
%     for i = 1 : length(BAH)
%         E_BAH(i) = F(BAH(i), ramp(r));
%     end
%     plot(BAH, E_BAH, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Bond Angle BAH");
% ylabel("F_{BAH}");
% fontsize(fs,"points");
% 
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% for r = 1 : length(ramp)
%     E_BAH = zeros(1,length(BAH));
%     for i = 1 : length(BAH)
%         E_BAH(i) = G(BAH(i), ramp(r));
%     end
%     plot(BAH, E_BAH, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Bond Angle BAH");
% ylabel("G_{BAH}");
% fontsize(fs,"points");
% 
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% for r = 1 : length(ramp)
%     Fs_AH = [0, 0.1, 3.2+(HA_max-3.3)*ramp(r), 3.3+(HA_max-3.3)*ramp(r)];
%     fs_AH = zeros(1,length(d_HA));
%     for i = 1 : length(d_HA)
%         fs_AH(i) = Fade(d_HA(i), Fs_AH);
%     end
%     plot(d_HA, fs_AH, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Distance HA");
% ylabel("Fade_{HA}");
% fontsize(fs,"points");
% 
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% AHD = linspace(0,pi,1000);
% for r = 1 : length(ramp)
%     F_AHD = [-0.95*ramp(r)-0.05, -0.95*ramp(r), 1, 1.05];
%     f_AHD = zeros(1,length(AHD));
%     for i = 1 : length(AHD)
%         f_AHD(i) = Fade(-cos(AHD(i)), F_AHD);
%     end
%     plot(AHD, f_AHD, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Bond Angle AHD");
% ylabel("Fade_{AHD}");
% fontsize(fs,"points");
% 
% figure;
% set(gcf, "color", "w");
% hold on;
% ramp = 0 : 0.2 : 1;
% BAH = linspace(0,pi,1000);
% for r = 1 : length(ramp)
%     F_BAH = [-0.4371*ramp(r)-0.562949, -0.4371*ramp(r), 1, 1.05];
%     f_BAH = zeros(1,length(BAH));
%     for i = 1 : length(AHD)
%         f_BAH(i) = Fade(-cos(BAH(i)), F_BAH);
%     end
%     plot(BAH, f_BAH, 'LineWidth', 2, "Color", [0.8,0.8,0.8]*ramp(r));
% end
% hold off;
% legend("ramp=0", "ramp=0.2", "ramp=0.4", "ramp=0.6", "ramp=0.8", "ramp=1");
% xlabel("Bond Angle BAH");
% ylabel("Fade_{BAH}");
% fontsize(fs,"points");
% 
% 
% figure;
% set(gcf, "color", "w");
% x = linspace(-3,3,1000);
% y = zeros(1,length(x));
% for i = 1 : length(x)
%     y(i) = f(x(i));
% end
% plot(x, y, 'LineWidth', 2, "Color", 'k');
% xlabel("x");
% ylabel("f");
% fontsize(fs,"points");
% 
% 
% function E_BAH = F(BAH, ramp)
% d = 0.75;
% m = 1.6;
% l = 0.357;
% 
% xmin = pi*(2/3-l);
% xmax = pi*2/3;
% xl = xmin*(1-ramp);
% 
% if BAH > xmax
%     F = d/2*cos(3*(pi-BAH)) + (d-1)/2;
% elseif BAH >= xl
%     z = xmin + (BAH-xl) / (xmax-xl) * (xmax-xmin);
%     F = m/2*cos(pi-(2/3*pi-z)/l) + (m-1)/2;
% else
%     F = m - 1/2;
% end
% 
% E_BAH = F;
% end
% 
% 
% function E_BAH = G(BAH, ramp)
% d = 0.75;
% m = 1.6;
% l = 0.357;
% 
% xmin = pi*(2/3-l);
% xmax = pi*2/3;
% xl = xmin*(1-ramp);
% 
% if BAH > xmax
%     G = d - 1/2;
% elseif BAH >= xl
%     z = xmin + (BAH-xl) / (xmax-xl) * (xmax-xmin);
%     G = (m-d)/2*cos(pi-(2/3*pi-z)/l) + (m+d-1)/2;
% else
%     G = m - 1/2;
% end
% 
% E_BAH = G;
% end



% hbond_bb: Backbone-backbone hbonds
% Donor is NH (hbdon_PBA), acceptor is CO (hbacc_PBA).
% Need also sequence separation to read HBEval.csv (headers in schema),
% then know which fade intervals from HBFadeIntervals to use, and which
% polynomial coefficients from HBPoly1D to use.
% hbonds_geom.cc from Rosetta implements this calculation.
% @ E_HA = F_AHD * F_BAH * poly(d_HA)
% @ E_AHD = Fs_AH * F_BAH * poly_s(AHD)
% @ E_B2BAH = E_B2BAH_calc(chi, BAH)
% @ WH strength of donor, WA strength of acceptor
% @ energy = Sum[WH*WA*f(E_HA+E_AHD+E_B2BAH)], f smooth min
% @ count: number of H-bonds
% @ oversat: number of oversaturated H-bonds, i.e. more than 2 hbonds to carbonyls
% @ coordinates: size 7nx3 matrix, 7n atoms (N, H, CA, 1HA, 2HA, C, O)
% @ ramp: this newly added ramp term is for relaxing the energy functions,
% so lenient H-bond requirements with wider energy basin regions when ramp
% is set high near 1, and the original strict H-bond requirements when ramp
% is set to 0.
function [energy, count, oversat] = E_hbond_ramping(coordinates, ramp)
Hbonds = [];

energy = 0;
count = 0; % number of hydrogen bonds
WH = 1.41;
WA = 1.08;
E_cutoff = -0.25;
threshold = 3.4473*2*sqrt(3);

HA_max = 6;

% Fade function parameters
Fs_AH = [0, 0.1, 3.2+(HA_max-3.3)*ramp, 3.3+(HA_max-3.3)*ramp];
F_AHD = [-0.95*ramp-0.05, -0.95*ramp, 1, 1.05];
F_BAH = [-0.4371*ramp-0.562949, -0.4371*ramp, 1, 1.05];

n = size(coordinates,1) / 7; % number of residues
for i = 1 : n % acceptor residue number
    for j = 1 : n % donor residue number
        if j ~= i
            % Compute atom-pair interaction for residues close enough
            if sum(abs(coordinates(7*i-4,:) - coordinates(7*j-4,:))) <= threshold

                % relevant atom coordinates
                B2 = coordinates(i*7-4,:); % Ca
                B = coordinates(i*7-1,:); % C
                A = coordinates(i*7,:); % O
                H = coordinates(j*7-5,:); % H
                D = coordinates(j*7-6,:); % N

                % calculate distance and angles
                d_HA = norm(A-H);
                AHD = acos(dot((A-H)/norm(A-H), (D-H)/norm(D-H)));
                BAH = acos(dot((B-A)/norm(B-A), (H-A)/norm(H-A)));
                chi = torsion(H, A, B, B2);

                % fade functions
                fs_AH = Fade(d_HA, Fs_AH);
                f_AHD = Fade(-cos(AHD), F_AHD);
                f_BAH = Fade(-cos(BAH), F_BAH);

                % polynomial parameters for E_HA, E_AHD
                xmin_HA = 1.384038127;
                xmax_HA = 2.998103943;
                fmin_HA = 1.1;
                fmax_HA = 1.1;
                C_HA = [-0.5307601, 6.47949946, -22.39522814, -55.14303544, 708.3094524, -2619.493182, 5227.88058, -6043.312116, 3806.046762, -1007.660241];
                xroot_HA = 2;

                xmin_AHD = 1.143564639;
                xmax_AHD = 3.1416;
                fmin_AHD = 1.1;
                fmax_AHD = 1.1;
                C_AHD = [0.47683259, -9.54524724, 83.62557693, -420.5586777, 1337.193549, -2786.262657, 3803.178227, -3278.628799, 1619.041162, -347.5015791];

                % calculate energies
                E_HA = f_AHD*f_BAH*polynomial_HA(d_HA, xmin_HA, xroot_HA, xmax_HA, HA_max, fmin_HA, fmax_HA, C_HA, ramp);
                E_AHD = fs_AH*f_BAH*polynomial_AHD(AHD, xmin_AHD, xmax_AHD, fmin_AHD, fmax_AHD, C_AHD, ramp);
                E_B2BAH = E_B2BAH_calc(chi, BAH, ramp);

                if fs_AH ~= 0 && f_AHD ~= 0 && f_BAH ~= 0
                    E = f(E_HA+E_AHD+E_B2BAH);
                    energy = energy + WH*WA*E;

                    if E < E_cutoff
                        fs_AH = Fade(d_HA, [0, 0.1, 3.2, 3.3]);
                        f_AHD = Fade(-cos(AHD), [-0.05, 0, 1, 1.05]);
                        f_BAH = Fade(-cos(BAH), [-0.562949, 0, 1, 1.05]);
                        if fs_AH ~= 0 && f_AHD ~= 0 && f_BAH ~= 0
                            E_HA = f_AHD*f_BAH*polynomial_HA(d_HA, xmin_HA, xroot_HA, xmax_HA, HA_max, fmin_HA, fmax_HA, C_HA, 0);
                            E_AHD = fs_AH*f_BAH*polynomial_AHD(AHD, xmin_AHD, xmax_AHD, fmin_AHD, fmax_AHD, C_AHD, 0);
                            E_B2BAH = E_B2BAH_calc(chi, BAH, 0);
                            E = f(E_HA+E_AHD+E_B2BAH);
                            if E < E_cutoff
                                count = count + 1;
                                Hbonds = [Hbonds; j, i];
                            end
                        end
                    end
                end
            end
        end
    end
end

[score, oversat] = Hbond_score(Hbonds, n);
energy = energy + score;
end


function [score, oversat] = Hbond_score(Hbonds, n)
score = 0;
oversat = 0;
penalty = 2;
reward = 2;

if Hbonds
    % assign rewards for long-range H-bonds
    hd = abs(Hbonds(:,1)-Hbonds(:,2));
    hd = min(hd, n-hd);
    score = score - reward*sum(log(max(hd-2, 1)));

    % % assign rewards for crossing H-bonds
    % cross = 0;
    % sortedHbonds = sort(Hbonds, 2);
    % sortedHbonds = sortrows(sortedHbonds, 1);
    % for i = 1 : size(sortedHbonds,1)
    %     for j = i+1 : size(sortedHbonds,1)
    %         if sortedHbonds(j,1) > sortedHbonds(i,1) && sortedHbonds(j,1) < sortedHbonds(i,2) && sortedHbonds(j,2) > sortedHbonds(i,2)
    %             cross = cross + 1;
    %         end
    %     end
    % end
    % score = score - sum(log(cross+1));

    % assign rewards for no consecutive residues without H-bonds
    sortedHbonds = sort(Hbonds(:));
    separation = [sortedHbonds(2:end); sortedHbonds(1)+n] - sortedHbonds;
    sigmoid = @(x) 1 ./ (1 + exp(x*1.5-6));
    score = score - reward*sigmoid(max(separation));

    % assign penalty for oversaturated Hbond acceptor, i.e. more than 2
    % hbonds to carbonyls
    acceptors = Hbonds(:,2);
    unique_elements = unique(acceptors);
    for u = 1 : length(unique_elements)
        oversat = oversat + max(sum(acceptors==unique_elements(u))-2, 0);
    end

    score = score + oversat*penalty;
end
end


function y = f(x)
if x < -0.1
    y = x;
elseif x < 0.1
    y = -0.025 + x/2 - 2.5*x^2;
else
    y = 0;
end
end


% Calculate energy E_B2BAH for SP2 hybrid
function energy = E_B2BAH_calc(chi, BAH, ramp)
d = 0.75;
m = 1.6;
l = 0.357;

H = (cos(2*chi)+1)/2;

xmin = pi*(2/3-l);
xmax = pi*2/3;
xl = xmin*(1-ramp);

if BAH > xmax
    F = d/2*cos(3*(pi-BAH)) + (d-1)/2;
elseif BAH >= xl
    z = xmin + (BAH-xl) / (xmax-xl) * (xmax-xmin);
    F = m/2*cos(pi-(2/3*pi-z)/l) + (m-1)/2;
else
    F = m - 1/2;
end

if BAH > xmax
    G = d - 1/2;
elseif BAH >= xl
    z = xmin + (BAH-xl) / (xmax-xl) * (xmax-xmin);
    G = (m-d)/2*cos(pi-(2/3*pi-z)/l) + (m+d-1)/2;
else
    G = m - 1/2;
end

energy = H*F + (1-H)*G;
end


% Fade function given I = [a,b,c,d]
% f=0 if x<=a or x>=d
% f=1 if b<=x<=c
% cubic Hermite spline function on [a,b] and [c,d] to ensure continuous derivative
function y = Fade(x, I)
a = I(1);
b = I(2);
c = I(3);
d = I(4);

if x <= a
    y = 0;
elseif x <= b
    z = (x-a) / (b-a);
    y = -2*z^3 + 3*z^2;
elseif x <= c
    y = 1;
elseif x <= d
    z = (x-c) / (d-c);
    y = 2*z^3 - 3*z^2 + 1;
else
    y = 0;
end
end


% Polynomial evaluation
% f=fmin if x<=xmin
% f=fmax if x>=xmax
% f=C[1]*x^(n-1) + C[2]*x^(n-2) + ... + C[n-1]*x + C[n] if xmin<x<xmax
function y = polynomial_HA(x, xmin, xroot, xmax, xb, fmin, fmax, C, ramp)
xl = xmin*(1-ramp);
xr = xmax+(xb-xmax)*ramp;

if x <= xl
    y = fmin;
elseif x <= xroot
    z = xmin + (x-xl) / (xroot-xl) * (xroot-xmin);
    y = C(1);
    for i = 2 : length(C)
        y = y*z + C(i);
    end
elseif x <= xr
    z = xroot + (x-xroot) / (xr-xroot) * (xmax-xroot);
    y = C(1);
    for i = 2 : length(C)
        y = y*z + C(i);
    end
else
    y = fmax;
end
end

% Polynomial evaluation
% f=fmin if x<=xmin
% f=fmax if x>=xmax
% f=C[1]*x^(n-1) + C[2]*x^(n-2) + ... + C[n-1]*x + C[n] if xmin<x<xmax
function y = polynomial_AHD(x, xmin, xmax, fmin, fmax, C, ramp)
xl = xmin*(1-ramp);

if x <= xl
    y = fmin;
elseif x <= xmax
    z = xmin + (x-xl) / (xmax-xl) * (xmax-xmin);
    y = C(1);
    for i = 2 : length(C)
        y = y*z + C(i);
    end
else
    y = fmax;
end
end

function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
m = cross(n1, b2/norm(b2));
x = dot(n1, n2);
y = dot(m, n2);
chi = atan2(y, x);
end
