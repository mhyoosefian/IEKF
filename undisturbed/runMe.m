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
Xh_iekf = cell(trials, 1);
Xh_ukf = cell(trials, 1);
X0 = zeros(3, 1);
m0 = 0.01*randn(3, 1);

n = size(X0,1);   % State vector size
m = 2;            % Measurement noise vector size
l = 3;            % Process noise vector size
L = l + m + n;

for iter=1:trials
    X{iter} = zeros(n, N);
    Xh_ekf{iter} = zeros(n ,N);
    Xh_iekf{iter} = zeros(n ,N);
    Xh_ukf{iter} = zeros(n ,N);
    X{iter}(:, 1) = X0;
    Xh_ekf{iter}(:, 1) = m0;
    Xh_iekf{iter}(:, 1) = m0;
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
    Q = 0.000001*eye(l);
    R = 2*eye(m);
    P = diag([.1 .1 (pi/50)]);
    P_ekf = P;
    P_ukf = P;
    W = [cos(Xh_iekf{iter}(3, t)) -sin(Xh_iekf{iter}(3, t)) 0; 
         sin(Xh_iekf{iter}(3, t)) cos(Xh_iekf{iter}(3, t)) 0;
         0 0 1];
    jacob_W_theta = [-sin(Xh_iekf{iter}(3, t)) -cos(Xh_iekf{iter}(3, t)) 0;
                     cos(Xh_iekf{iter}(3, t)) -sin(Xh_iekf{iter}(3, t)) 0;
                     0 0 0];
    P_iekf = W'*P*W + P(3, 3)*jacob_W_theta'*P*jacob_W_theta;
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
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% IEKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = A_k(u, dt);          % jacobian of 'f' wrt 'X'
        G = eye(3);          % jacobian of 'f' wrt 'w'
        
        % 1- Propagation
        Xh_minus_iekf = f(Xh_iekf{iter}(:, t), u, zeros(l,1), dt);
        T = [cos(Xh_minus_iekf(3)) sin(Xh_minus_iekf(3)); -sin(Xh_minus_iekf(3)) cos(Xh_minus_iekf(3))];
        jacob_T_theta = [-sin(Xh_minus_iekf(3)) cos(Xh_minus_iekf(3)); -cos(Xh_minus_iekf(3)) -sin(Xh_minus_iekf(3))];
        W = [cos(Xh_minus_iekf(3)) -sin(Xh_minus_iekf(3)) 0; 
             sin(Xh_minus_iekf(3)) cos(Xh_minus_iekf(3)) 0;
             0 0 1];
        jacob_W_theta = [-sin(Xh_minus_iekf(3)) -cos(Xh_minus_iekf(3)) 0;
                         cos(Xh_minus_iekf(3)) -sin(Xh_minus_iekf(3)) 0;
                         0 0 0];

        Q_iekf = W'*Q*W + P_iekf(3, 3)*jacob_W_theta'*Q*jacob_W_theta;
        P_minus_iekf = F*P_iekf*F.' + G*Q_iekf*G.';
        
        H = [1 0 0; 0 1 0];     % jacobian of 'h' wrt 'X'
        M = eye(m);                   % jacobian of 'h' wrt 'v'
        
        % 2- Update
        R_iekf = T*R*T' + P_minus_iekf(3, 3)*jacob_T_theta*R*jacob_T_theta';
        K = P_minus_iekf*H.'*inv(H*P_minus_iekf*H.' + M*R_iekf*M.');
        Xh_iekf{iter}(:,t+1) = Xh_minus_iekf + W*K*T*(y - h(Xh_minus_iekf, zeros(m,1)));
        P_iekf = (eye(n) - K*H)*P_minus_iekf;


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
rmse_iekf = zeros(n, N);
rmse_ukf = zeros(n, N);
for i=1:N
    for j=1:n
        tmp_ekf = 0;
        tmp_iekf = 0;
        tmp_ukf = 0;
        for k=1:trials
            tmp_ekf = tmp_ekf + (X{k}(j, i) - Xh_ekf{k}(j, i))^2;
            tmp_iekf = tmp_iekf + (X{k}(j, i) - Xh_iekf{k}(j, i))^2;
            tmp_ukf = tmp_ukf + (X{k}(j, i) - Xh_ukf{k}(j, i))^2;
        end
        rmse_ekf(j, i) = sqrt(tmp_ekf/trials);
        rmse_iekf(j, i) = sqrt(tmp_iekf/trials);
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
plot(dt:dt:N*dt, rmse_iekf(1,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(1,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $x$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(2,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf(2,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(2,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $y$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex')
grid on
box on

figure;
plot(dt:dt:N*dt, rmse_ekf(3,:), 'b-.', 'LineWidth', 1);
hold on
plot(dt:dt:N*dt, rmse_iekf(3,:), 'r-.', 'LineWidth', 1);
plot(dt:dt:N*dt, rmse_ukf(3,:), 'm-.', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('RMSE of $\theta$', 'FontSize', 17, 'Interpreter', 'latex');
legend('EKF', 'IEKF', 'UKF', 'FontSize', 17, 'Interpreter', 'latex');
set(gca,'FontSize', 14, 'TickLabelInterpreter', 'latex')
grid on
box on
