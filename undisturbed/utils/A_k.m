function out = A_k(u, dt)

v = u(1);
omega = u(2);
out = [0 omega 0;
       -omega 0 v;
       0 0 0];
