% Function to find a feasible direction
function [p_f, angles_f] = feasible(p, angles, Boundaries)
p_f = p;
angles_f = angles + p_f;
angles_f = mod(angles_f+pi, 2*pi) - pi;
in_or_not = double(inpolygon(angles_f(1:2:end-1), angles_f(2:2:end), Boundaries(1,:), Boundaries(2,:)).');
in_or_not = repelem(in_or_not, 2);

while sum(in_or_not==0)>0 && sum(abs(p_f(in_or_not==0)))>0.001
    in_or_not(in_or_not==0) = -0.8;
    p_f = p_f.*in_or_not.';
    angles_f = angles + p_f;
    angles_f = mod(angles_f+pi, 2*pi) - pi;
    in_or_not = double(inpolygon(angles_f(1:2:end-1), angles_f(2:2:end), Boundaries(1,:), Boundaries(2,:)).');
    in_or_not = repelem(in_or_not, 2);
end
end