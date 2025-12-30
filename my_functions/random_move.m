% Function to generate random move direction
function p = random_move(k, n)

a = 2*pi*rand(1,n);
r = k*sqrt(rand(1,n));
p = reshape([r.*sin(a); r.*cos(a)], [2*n,1]);

end