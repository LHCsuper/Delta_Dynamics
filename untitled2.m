clc
clear
close all


%% 各个参数
mm = 1e-3;
I1 =  0.065490;
I_motor = 1.303;%0.002085
m1 = 1.672;
m2 = 0.11;
mp = 0.274;

%% 仿真

%轨迹生成
%% ===== 定义分拣点 =====
A = [-0.3, 0, -0.7];
B = [-0.3, 0, -0.65];
C = [-0.2, 0, -0.6];
D = [ 0.0, 0, -0.6];
E = [ 0.2, 0, -0.6];
F = [ 0.3, 0, -0.65];
G = [ 0.3, 0, -0.7];

T  = 1.5 ;
dt = 0.001;
t  = 0:dt:T;
nCycle=1;
%% ===== 生成轨迹 =====
traj = generateSortingTrajectory(A,B,C,D,E,F,G,T,t,nCycle);
fprintf('theta start = [% .6f % .6f % .6f]\n', traj.theta(:,1));
fprintf('theta end   = [% .6f % .6f % .6f]\n', traj.theta(:,end));

fprintf('P start = [% .6f % .6f % .6f]\n', traj.x(1), traj.y(1), traj.z(1));
fprintf('P end   = [% .6f % .6f % .6f]\n', traj.x(end), traj.y(end), traj.z(end));
%% ===== 数值计算 v,a =====
vx = gradient(traj.x, dt);
vy = gradient(traj.y, dt);
vz = gradient(traj.z, dt);

ax = gradient(vx, dt);
ay = gradient(vy, dt);
az = gradient(vz, dt);

%% ===== 逆运动学 =====
N = length(t);
theta = zeros(3,N);

for i = 1:N
    P = [traj.x(i);
         traj.y(i);
         traj.z(i)];
    theta(:,i) = deltaIK(P);
end

%% ===== 关节速度加速度 =====
dtheta  = gradient(theta, dt);
ddtheta = gradient(dtheta, dt);

%% 控制反证
%% ================= 参数范围 =================
kp_list = 20:10:120;     % 位置增益扫描范围
kd_list = 10:5:60;       % 速度增益扫描范围

n_kp = length(kp_list);
n_kd = length(kd_list);

RMSE_map = zeros(n_kd, n_kp);  % 存储结果
%% ================= 控制参数 =================
for i = 1:n_kd
    for j = 1:n_kp
        
        kp = kp_list(j);
        kd = kd_list(i);

KP = diag([kp,kp,kp]);     % 位置增益
KD = diag([kd,kd,kd]);       % 速度增益

n = length(traj.t);
step = traj.t(2) - traj.t(1);

%% ================= 初始化 =================
q   = zeros(3,n);
dq  = zeros(3,n);
tau = zeros(3,n);

% 初始状态与轨迹一致
q(:,1)  = traj.theta(:,1);
dq(:,1) = traj.dtheta(:,1);

%% ================= 主循环 =================
for k = 1:n-1
    
    %% ---------- 期望 ----------
    qd   = traj.theta(:,k);
    dqd  = traj.dtheta(:,k);
    ddqd = traj.ddtheta(:,k);
    
    %% ---------- 当前 ----------
    qk  = q(:,k);
    dqk = dq(:,k);
    
    %% ---------- 误差 ----------
    e  = qd  - qk;
    de = dqd - dqk;
    
    %% ========== 1️⃣ 当前几何一致的末端位置 ==========
    Pk = deltaFK(qk);    % 必须使用当前 qk
    xk=Pk(1);
    yk=Pk(2);
    zk=Pk(3);
    %% ========== 2️⃣ 生成控制加速度 ==========
    ddq_cmd = ddqd + KD*de + KP*e;
    
    %% ========== 3️⃣ 逆动力学得到控制力矩 ==========
    %计算力矩法
    tau(:,k) = tau_delta( ...
        I1, I_motor, ...
        0, 0, -9.81, ...          % 重力方向
        ddq_cmd(1), ddq_cmd(2), ddq_cmd(3), ...
        dqk(1), dqk(2), dqk(3), ...
        m1, m2, mp, ...
        qk(1), qk(2), qk(3), ...
        xk, yk, zk);
   
    %纯PD法
        % tau(:,k) = KP*e + KD*de;
    %% 正动力学更新位置
    ddq = dynamics_inverse( ...
            tau(:,k), ...
            qk, dqk, ...
            I1, I_motor, ...
            0, 0, -9.81, ...
            m1, m2, mp, ...
            xk, yk, zk);
    
    %% ========== 5️⃣ 半隐式 Euler 积分 ==========
    dq(:,k+1) = dqk + ddq * step;
    q(:,k+1)  = qk  + dq(:,k+1) * step;
    
end


t_plot = traj.t;

theta_plot  = traj.theta;
dtheta_plot = traj.dtheta;
ddtheta_plot = traj.ddtheta;



% for i = 1:3
% 
%     figure;   % 每次循环新建一个窗口
% 
%     % 先画 Actual（蓝色实线）
%     plot(traj.t, q(i,:), 'b', 'LineWidth', 1.2); 
%     hold on;
% 
%     % 后画 Desired（红色虚线）
%     plot(traj.t, traj.theta(i,:), 'r--', 'LineWidth', 1.8);
% 
%     grid on;
%     xlabel('Time (s)');
%     ylabel(['\theta_', num2str(i)]);
%     legend('MATLAB','Adams&Simulink');
%     title(['Joint ', num2str(i), ' Tracking Performance']);
% 
% end

% ===== 计算实际末端位置 =====
N = length(traj.t);

x_actual = zeros(1,N);
y_actual = zeros(1,N);
z_actual = zeros(1,N);

for k = 1:N
    P = deltaFK(q(:,k));   % 正运动学
    x_actual(k) = P(1);
    y_actual(k) = P(2);
    z_actual(k) = P(3);
end

% figure;
% plot3(traj.x, traj.y, traj.z, 'r--', 'LineWidth', 2); hold on;
% plot3(x_actual, y_actual, z_actual, 'b', 'LineWidth', 1.5);
% 
% grid on;
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Z (m)');
% legend('Desired','Actual');
% title('Cartesian Space Tracking');
% axis equal;
% 计算瞬时欧氏误差
error_vec = sqrt( ...
    (traj.x - x_actual).^2 + ...
    (traj.y - y_actual).^2 + ...
    (traj.z - z_actual).^2 );

% 平均误差
mean_error = mean(error_vec);

% 均方根误差（更常用）
rmse_error = sqrt(mean(error_vec.^2));
   RMSE_map(i,j) = rmse_error;
   fprintf('kp=%.1f, kd=%.1f, RMSE=%.6f\n', kp, kd, rmse_error);
    end
end

% 展开为列向量
[Kp_grid, Kd_grid] = meshgrid(kp_list, kd_list);

Kp_vec = Kp_grid(:);
Kd_vec = Kd_grid(:);
E_vec  = RMSE_map(:);

% 构造设计矩阵
X = [ ...
    ones(size(Kp_vec)), ...
    Kp_vec, ...
    Kd_vec, ...
    Kp_vec.^2, ...
    Kd_vec.^2, ...
    Kp_vec .* Kd_vec];

% 最小二乘求解
beta = X \ E_vec;

% 系数
b0 = beta(1);
b1 = beta(2);
b2 = beta(3);
b3 = beta(4);
b4 = beta(5);
b5 = beta(6);

fprintf('Fitted model:\n');
fprintf('E = %.6e + %.6e*Kp + %.6e*Kd + %.6e*Kp^2 + %.6e*Kd^2 + %.6e*Kp*Kd\n', ...
    b0,b1,b2,b3,b4,b5);

figure;
imagesc(kp_list, kd_list, RMSE_map);
set(gca,'YDir','normal');
xlabel('K_p');
ylabel('K_d');
title('RMSE Heatmap');
colorbar;


% 最大误差
max_error = max(error_vec);

fprintf('平均跟踪误差 = %.6f m\n', mean_error);
fprintf('RMSE误差      = %.6f m\n', rmse_error);
fprintf('最大跟踪误差 = %.6f m\n', max_error);

figure;
clf;
hold on;

% θ1 — 实线
plot(t_plot, theta_plot(1,:), 'b-',  'LineWidth', 2);
% θ3 — 实线
plot(t_plot, theta_plot(3,:), 'k-',  'LineWidth', 2);
% θ2 — 虚线（避免遮挡）
plot(t_plot, theta_plot(2,:), 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Joint Angle (rad)');
title('Joint Angular Displacements');

legend({'$\theta_1$', ...
        '$\theta_2$', ...
        '$\theta_3$'}, ...
        'Interpreter','latex');





% figure;
% clf;
% hold on;
% % θ1 — 实线
% plot(t_plot, dtheta_plot(1,:), 'b-',  'LineWidth', 2);
% 
% 
% 
% % θ3 — 实线
% plot(t_plot, dtheta_plot(3,:), 'k-',  'LineWidth', 2);
% % θ2 — 虚线（解决遮挡）
% plot(t_plot, dtheta_plot(2,:), 'r--', 'LineWidth', 2);
% grid on;
% xlabel('Time (s)');
% ylabel('Angular Velocity (rad/s)');
% title('Joint Angular Velocities');
% 
% legend({'$\dot{\theta}_1$', ...
%         '$\dot{\theta}_2$', ...
%         '$\dot{\theta}_3$'}, ...
%         'Interpreter','latex');

figure;
clf;
hold on;
% θ1 — 实线
plot(t_plot, ddtheta_plot(1,:), 'b-',  'LineWidth', 2);

% θ3 — 实线
plot(t_plot, ddtheta_plot(3,:), 'k-',  'LineWidth', 2);
% θ2 — 虚线（避免遮挡）
plot(t_plot, ddtheta_plot(2,:), 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Angular Acceleration (rad/s^2)');
title('Joint Angular Accelerations');

legend({'$\ddot{\theta}_1$', ...
        '$\ddot{\theta}_2$', ...
        '$\ddot{\theta}_3$'}, ...
        'Interpreter','latex');



