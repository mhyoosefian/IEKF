function out1 = f(X, u, w, dt)
% INPUT: 1- X (n by 1): state vector
%        2- u (m by 1): controller
%        3- w (l by 1): process noise

% OUTPUT: out (n by 1): X(k+1) = f(X(k), u(k), w(k))


out = zeros(length(X), 1);
v = u(1);
omega = u(2);
theta = X(3);
out(1) = v*cos(theta);
out(2) = v*sin(theta);
out(3) = omega;
out = out + w;

out1 = dt*out + X;
