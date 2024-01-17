function y = h(X, v)
% INPUT: 1- X (n by 1): state vector
%        2- v (m by 1): process noise

% OUTPUT: out (n by 1): y(k) = h(X(k), v(k))

% n = size(X);
% m = size(v);
% 
% N1 = 20;
% E1 = 0;
% N2 = 0;
% E2 = 20;
% 
% 
% y = [sqrt((X(1)-N1)^2 + (X(2)-E1)^2); ...
% sqrt((X(1)-N2)^2 + (X(2)-E2)^2)] + v;

y = [eye(2) zeros(2, 5)]*X + v;