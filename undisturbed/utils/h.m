function y = h(X, v)
% INPUT: 1- X (n by 1): state vector
%        2- v (m by 1): process noise

% OUTPUT: out (n by 1): y(k) = h(X(k), v(k))

y = [1 0 0; 0 1 0]*X + v;
% y = X + v;