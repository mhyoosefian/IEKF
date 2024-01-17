function out = A_k1(X, u, dt)

v = u(1);
omega = u(2);
theta = X(3);
Cx = [1 0 0 1];
Cy = [0 -1 1 0];
A = [0 0 1 1;
    0 0 0 0;
    -1 -1 0 0;
    0 0 0 0];
out = [0 omega 0 Cx*cos(theta)+Cy*sin(theta);
       -omega 0 v -Cx*sin(theta)+Cy*cos(theta);
       0 0 0 zeros(1, 4)
       zeros(4, 1) zeros(4, 1) zeros(4, 1) A];
