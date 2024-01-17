function F = jacob_f_X(X,u, dt)
v = u(1);
theta = X(3);
Cx = [1 0 0 1];
Cy = [0 -1 1 0];
A = [0 0 1 1;
    0 0 0 0;
    -1 -1 0 0;
    0 0 0 0];
F = eye(7) + dt*[0 0 -v*sin(theta) Cx;
    0 0 v*cos(theta) Cy;
    0 0 0 zeros(1, 4);
    zeros(4, 1) zeros(4, 1) zeros(4, 1) A];