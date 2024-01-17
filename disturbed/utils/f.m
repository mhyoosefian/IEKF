function out = f(X, u, w, dt)
% INPUT: 1- X (n by 1): state vector
%        2- u (m by 1): controller
%        3- w (l by 1): process noise

% OUTPUT: out (n by 1): X(k+1) = f(X(k), u(k), w(k))

% n = size(X);
% m = size(u);
% l = size(w);
% 
% T = 0.1;
% Phi = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
% B = [1 1 -1 -1].';
% out = Phi*X + B*u + w;

out = zeros(length(X), 1);
Cx = [1 0 0 1];
Cy = [0 -1 1 0];
v = u(1);
omega = u(2);
A = [0 0 1 1;
    0 0 0 0;
    -1 -1 0 0;
    0 0 0 0];
theta = X(3);
d = X(4:end);
out(1) = v*cos(theta) + Cx*d;
out(2) = v*sin(theta) + Cy*d;
out(3) = omega;
out(4:end) = A*d;
out = out + w;

out = dt*out + X;
