function out = A_k2(u, dt)

v = u(1);
omega = u(2);
Cx = [1 0 0 1];
Cy = [0 -1 1 0];
A = [0 0 1 1;
    0 0 0 0;
    -1 -1 0 0;
    0 0 0 0];
temp = [0 omega; -omega 0];

Aw = A + kron(temp, eye(2));

out = [0 omega 0 Cx;
       -omega 0 v Cy;
       0 0 0 zeros(1, 4)
       zeros(4, 1) zeros(4, 1) zeros(4, 1) Aw];
