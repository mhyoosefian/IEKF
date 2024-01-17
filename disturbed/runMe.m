% This is the code that you should run :)


%%

clc;
clear all
addpath('utils')
dt = 0.05;
N = 200/dt;
trials = 10;
rng(1);

X = cell(trials, 1);
Xh_ekf = cell(trials, 1);
Xh_iekf1 = cell(trials, 1);
Xh_iekf2 = cell(trials, 1);
Xh_ukf = cell(trials, 1);
X0 = zeros(7, 1);
m0 = 0.0001*randn(7, 1);

n = size(X0,1);   % State vector size
m = 2;            % Measurement noise vector size
l = 7;            % Process noise vector size
L = l + m + n;

for iter=1:trials
    X{iter} = zeros(n, N);
    Xh_ekf{iter} = zeros(n ,N);
    Xh_iekf1{iter} = zeros(n ,N);
    Xh_iekf2{iter} = zeros(n ,N);
    Xh_ukf{iter} = zeros(n ,N);
    X{iter}(:, 1) = X0;
    Xh_ekf{iter}(:, 1) = m0;
    Xh_iekf1{iter}(:, 1) = m0;
    Xh_iekf2{iter}(:, 1) = m0;
    Xh_ukf{iter}(:, 1) = m0;
end
alpha = 1;          % Primary scaling parameter
beta = 2;           % Secondary scaling parameter
kappa = 0;          % Teritary scaling parameter
lambda = alpha^2*(L+kappa) - L;
etha_m = [lambda/(L+lambda); ones(2*L,1)*1/(2*(L+lambda))];
etha_c = [lambda/(L+lambda) + 1 - alpha^2 + beta; ones(2*L,1)*1/(2*(L+lambda))];

%%
for iter=1:trials
    t = 1;
    disp(iter);
    Q = 1e-5*eye(l);
    R = 2*eye(m);
    P = diag([.1 .1 (pi/20) .1 .1 .1 .1]);
    P_ekf = P;
    P_ukf = P;

    T = [cos(Xh_iekf1{iter}(3, t)) sin(Xh_iekf1{iter}(3, t)); -sin(Xh_iekf1{iter}(3, t)) cos(Xh_iekf1{iter}(3, t))];
    W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
        zeros(4, 3) eye(4)];
    tmp = [-sin(Xh_iekf1{iter}(3, t)) -cos(Xh_iekf1{iter}(3, t)); cos(Xh_iekf1{iter}(3, t)) -sin(Xh_iekf1{iter}(3, t))];
    jacob_W_theta = [[tmp zeros(2, 1); zeros(1, 3)] zeros(3, 4);
                     zeros(4, 3) zeros(4, 4)];
    P_iekf1 = W'*P*W + jacob_W_theta'*P(3, 3)*jacob_W_theta;
    
    T = [cos(Xh_iekf2{iter}(3, t)) sin(Xh_iekf2{iter}(3, t)); -sin(Xh_iekf2{iter}(3, t)) cos(Xh_iekf2{iter}(3, t))];
    W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
          zeros(4, 3) kron(T', eye(2))];
    tmp = [-sin(Xh_iekf2{iter}(3, t)) -cos(Xh_iekf2{iter}(3, t)); cos(Xh_iekf2{iter}(3, t)) -sin(Xh_iekf2{iter}(3, t))];
    jacob_W_theta = [[tmp zeros(2, 1); zeros(1, 3)] zeros(3, 4);
                     zeros(4, 3) kron(tmp, eye(2))];
    P_iekf2 = W'*P*W + jacob_W_theta'*P(3, 3)*jacob_W_theta;
    w = sqrt(Q)*randn(l, N);
    v = sqrt(R)*randn(m, N);
    [~,sqrt_Q,~] = ldl(Q);
    sqrt_R = chol(R, 'lower');
    
    while t<N
        u = [1*sin(2*pi*t/N); 4*2*pi/180];     % System input (u = [1.3*sin(2*pi*t/N); ...])
        
        X{iter}(:,t+1) = f(X{iter}(:,t), u, w(:,t), dt);     % True states
        y = h(X{iter}(:,t), v(:,t));                      % Observation
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = jacob_f_X(Xh_ekf{iter}(:,t), u, dt);          % jacobian of 'f' wrt 'X'
        G = jacob_f_w(Xh_ekf{iter}(:,t), u, dt);          % jacobian of 'f' wrt 'w'
        
        % 1- Propagation
        Xh_minus_ekf = f(Xh_ekf{iter}(:,t), u, zeros(l,1), dt);
        P_minus_ekf = F*P_ekf*F.' + G*Q*G.';
        
        H = jacob_h_X(Xh_minus_ekf);     % jacobian of 'h' wrt 'X'
        M = jacob_h_v(Xh_minus_ekf);     % jacobian of 'h' wrt 'v'
        
        % 2- Update
        K = P_minus_ekf*H.'*inv(H*P_minus_ekf*H.' + M*R*M.');
        Xh_ekf{iter}(:,t+1) = Xh_minus_ekf + K*(y - h(Xh_minus_ekf, zeros(m,1)));
        P_ekf = (eye(n) - K*H)*P_minus_ekf;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% IEKF1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = A_k1(Xh_iekf1{iter}(:,t), u, dt);          % jacobian of 'f' wrt 'X'
        G = eye(l);          % jacobian of 'f' wrt 'w'
        
        % 1- Propagation
        Xh_minus_iekf1 = f(Xh_iekf1{iter}(:,t), u, zeros(l,1), dt);
        T = [cos(Xh_minus_iekf1(3)) sin(Xh_minus_iekf1(3)); -sin(Xh_minus_iekf1(3)) cos(Xh_minus_iekf1(3))];
        W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
            zeros(4, 3) eye(4)];
        tmp = [-sin(Xh_minus_iekf1(3)) -cos(Xh_minus_iekf1(3)); cos(Xh_minus_iekf1(3)) -sin(Xh_minus_iekf1(3))];
        jacob_W_theta = [[tmp zeros(2, 1); zeros(1, 3)] zeros(3, 4);
                         zeros(4, 3) zeros(4, 4)];
        Q_iekf1 = W'*Q*W + P_iekf1(3, 3)*jacob_W_theta'*Q*jacob_W_theta;
        P_minus_iekf1 = F*P_iekf1*F.' + G*Q_iekf1*G.';
        
        H = [eye(2) zeros(2, 5)];     % jacobian of 'h' wrt 'X'
        M = eye(m);                   % jacobian of 'h' wrt 'v'
        
        % 2- Update
        T = [cos(Xh_minus_iekf1(3)) sin(Xh_minus_iekf1(3)); -sin(Xh_minus_iekf1(3)) cos(Xh_minus_iekf1(3))];
        W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
            zeros(4, 3) eye(4)];
        jacob_T_theta = [-sin(Xh_minus_iekf1(3)) cos(Xh_minus_iekf1(3)); -cos(Xh_minus_iekf1(3)) -sin(Xh_minus_iekf1(3))];
        R_iekf1 = T'*R*T + P_minus_iekf1(3, 3)*jacob_T_theta*R*jacob_T_theta';
        K = P_minus_iekf1*H.'*inv(H*P_minus_iekf1*H.' + M*R_iekf1*M.');
        Xh_iekf1{iter}(:,t+1) = Xh_minus_iekf1 + W*K*T*(y - h(Xh_minus_iekf1, zeros(m,1)));
        P_iekf1 = (eye(n) - K*H)*P_minus_iekf1;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% IEKF2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = A_k2(u, dt);          % jacobian of 'f' wrt 'X'
        G = eye(l);          % jacobian of 'f' wrt 'w'
        
        % 1- Propagation
        Xh_minus_iekf2 = f(Xh_iekf2{iter}(:,t), u, zeros(l,1), dt);
        T = [cos(Xh_minus_iekf2(3)) sin(Xh_minus_iekf2(3)); -sin(Xh_minus_iekf2(3)) cos(Xh_minus_iekf2(3))];
        W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
            zeros(4, 3) kron(T.', eye(2))];
        tmp = [-sin(Xh_minus_iekf2(3)) -cos(Xh_minus_iekf2(3)); cos(Xh_minus_iekf2(3)) -sin(Xh_minus_iekf2(3))];
        jacob_W_theta = [[tmp zeros(2, 1); zeros(1, 3)] zeros(3, 4);
                         zeros(4, 3) kron(tmp, eye(2))];
        Q_iekf2 = W'*Q*W + P_iekf2(3, 3)*jacob_W_theta'*Q*jacob_W_theta;
        P_minus_iekf2 = F*P_iekf2*F.' + G*Q_iekf2*G.';
        
        H = [eye(2) zeros(2, 5)];     % jacobian of 'h' wrt 'X'
        M = eye(m);                   % jacobian of 'h' wrt 'v'
        
        % 2- Update
        T = [cos(Xh_minus_iekf2(3)) sin(Xh_minus_iekf2(3)); -sin(Xh_minus_iekf2(3)) cos(Xh_minus_iekf2(3))];
        W = [[T.' zeros(2, 1); zeros(1, 2) 1] zeros(3, 4);
            zeros(4, 3) kron(T.', eye(2))];
        jacob_T_theta = [-sin(Xh_minus_iekf2(3)) cos(Xh_minus_iekf2(3)); -cos(Xh_minus_iekf2(3)) -sin(Xh_minus_iekf2(3))];
        R_iekf2 = T'*R*T + P_minus_iekf2(3, 3)*jacob_T_theta*R*jacob_T_theta';
        K = P_minus_iekf2*H.'*inv(H*P_minus_iekf2*H.' + M*R_iekf2*M.');
        Xh_iekf2{iter}(:,t+1) = Xh_minus_iekf2 + W*K*T*(y - h(Xh_minus_iekf2, zeros(m,1)));
        P_iekf2 = (eye(n) - K*H)*P_minus_iekf2;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xha = [Xh_ukf{iter}(:,t); zeros(l,1); zeros(m,1)];     % Augmented estimation
        sqrt_P = chol(P_ukf,'lower');
        sqrt_Pa = [sqrt_P zeros(n,l) zeros(n,m)
                   zeros(l,n) sqrt_Q zeros(l,m)
                   zeros(m,n) zeros(m,l) sqrt_R];
        chi1 = [Xha, Xha*ones(1,L) + sqrt(L+lambda)*sqrt_Pa, Xha*ones(1,L) - sqrt(L+lambda)*sqrt_Pa];     % Sigma points
        chi1_X = chi1(1:n,:);
        chi1_W = chi1(n+1:n+l,:);
        chi1_V = chi1(n+l+1:L,:);
        for i=1:2*L+1
            chi2_X(:,i) = f(chi1_X(:,i), u, chi1_W(:,i), dt);     % First unsceted transform
        end
        
        % 1- Propagation
        Xh_minus_ukf = chi2_X*etha_m;
        P_minus_ukf = zeros(n,n);
        for i =1:2*L+1
            P_minus_ukf = P_minus_ukf + etha_c(i)*((chi2_X(:,i) - Xh_minus_ukf)*(chi2_X(:,i) - Xh_minus_ukf).');
        end
        Psi = zeros(m, 2*L+1);
        for i=1:2*L+1
            Psi(:,i) = h(chi2_X(:,i), chi1_V(:,i));     % Second unscented transform
        end
        yhat = Psi*etha_m;
        Pyy = zeros(m,m);
        Pxy = zeros(n,m);
        for i =1:2*L+1
            Pyy = Pyy + etha_c(i)*(Psi(:,i) - yhat)*(Psi(:,i) - yhat).';
            Pxy = Pxy + etha_c(i)*(chi2_X(:,i) - Xh_minus_ukf)*(Psi(:,i) - yhat).';
        end
        % 2- Update
        K = Pxy/Pyy;
        Xh_ukf{iter}(:,t+1) = Xh_minus_ukf + K*(y - yhat);
        P_ukf = P_minus_ukf - K*Pyy*K.';
        t = t+1;
    end
end

%%
rmse_ekf = zeros(n, N);
rmse_iekf1 = zeros(n, N);
rmse_iekf2 = zeros(n, N);
rmse_ukf = zeros(n, N);
for i=1:N
    for j=1:n
        tmp_ekf = 0;
        tmp_iekf1 = 0;
        tmp_iekf2 = 0;
        tmp_ukf = 0;
        for k=1:trials
            tmp_ekf = tmp_ekf + (X{k}(j, i) - Xh_ekf{k}(j, i))^2;
            tmp_iekf1 = tmp_iekf1 + (X{k}(j, i) - Xh_iekf1{k}(j, i))^2;
            tmp_iekf2 = tmp_iekf2 + (X{k}(j, i) - Xh_iekf2{k}(j, i))^2;
            tmp_ukf = tmp_ukf + (X{k}(j, i) - Xh_ukf{k}(j, i))^2;
        end
        rmse_ekf(j, i) = sqrt(tmp_ekf/trials);
        rmse_iekf1(j, i) = sqrt(tmp_iekf1/trials);
        rmse_iekf2(j, i) = sqrt(tmp_iekf2/trials);
        rmse_ukf(j, i) = sqrt(tmp_ukf/trials);
    end
end

%%
% figure;
% plot(dt:dt:N*dt, X{1}(1,:), 'k--', 'LineWidth', 1);
% hold on
% plot(dt:dt:N*dt, Xh_ekf{1}(1,:), 'b-.', 'LineWidth', 1);
% plot(dt:dt:N*dt, Xh_iekf1{1}(1,:), 'r-.', 'LineWidth', 1);
% plot(dt:dt:N*dt, Xh_ukf{1}(1,:), 'm-.', 'LineWidth', 1);
% xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
% ylabel('x (m)', 'FontSize', 17, 'Interpreter', 'latex');
% legend('Ground truth', 'EKF', 'IEKF1', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
% grid on
% box on

figure;
plot(dt:dt:N*dt, rmse_ekf(1,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(1,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(1,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(1,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $x$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(2,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(2,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(2,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(2,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $y$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(3,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(3,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(3,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(3,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $\theta$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(4,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(4,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(4,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(4,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $d_1$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(5,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(5,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(5,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(5,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $d_2$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(6,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(6,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(6,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(6,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $d_3$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(7,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf1(7,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_iekf2(7,:), 'k-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(7,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $d_4$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF1', 'IEKF2', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on
