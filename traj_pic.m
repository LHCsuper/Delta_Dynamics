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

T  = 0.25;
nTest = 125;

T_list    = zeros(1,nTest);   % 存储每次的T
rmse_list = zeros(1,nTest);   % 存储每次的RMSE
for ii = 1:nTest
T=T+0.01;
fprintf('周期T      = %.6f s\n', T);
T_list(ii) = T;
dt = 0.001;
t  = 0:dt:T;
nCycle=2;
%% ===== 生成轨迹 =====
traj = generateSortingTrajectory(A,B,C,D,E,F,G,T,t,nCycle);
% fprintf('theta start = [% .6f % .6f % .6f]\n', traj.theta(:,1));
% fprintf('theta end   = [% .6f % .6f % .6f]\n', traj.theta(:,end));
% 
% fprintf('P start = [% .6f % .6f % .6f]\n', traj.x(1), traj.y(1), traj.z(1));
% fprintf('P end   = [% .6f % .6f % .6f]\n', traj.x(end), traj.y(end), traj.z(end));
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
%% ================= 控制参数 =================
kp=100;
kd=60;

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

% 均方根误差（更常用）
rmse_error = sqrt(mean(error_vec.^2));

fprintf('RMSE误差      = %.6f m\n', rmse_error);

    rmse_list(ii) = rmse_error;   % 存储
end
% figure;
% plot(T_list, rmse_list, 'o-', 'LineWidth', 2);
% grid on;
% xlabel('Period T (s)');
% ylabel('RMSE (m)');
% title('Effect of Period on Tracking Accuracy');
% 取对数
X = log(T_list(:));
Y = log(rmse_list(:));

% 线性回归
p = polyfit(X, Y, 1);

c1 = p(1);   % 斜率
c0 = p(2);   % 截距

b = -c1;
a = exp(c0);

fprintf('Fitted model: E(T) = %.6e * T^(-%.4f)\n', a, b);

% 生成拟合曲线
T_fit = linspace(min(T_list), max(T_list), 200);
E_fit = a * T_fit.^(-b);

% 画图对比
figure('Color','w');
plot(T_list, rmse_list, 'o-', ...
    'LineWidth',2, ...
    'Color',[0 0.4470 0.7410]); 
hold on;

grid on;

% 坐标轴统一
set(gca, ...
    'FontSize',18, ...
    'LineWidth',1.5, ...
    'TickDir','in', ...
    'Box','on');

xlabel('T (s)', ...
    'FontSize',20, ...
    'FontWeight','bold');

ylabel('RMSE (m)', ...
    'FontSize',20, ...
    'FontWeight','bold');

title('T–RMSE Relationship', ...
    'FontSize',21);

legend('Data', ...
    'FontSize',16, ...
    'Location','best');

figure('Color','w');

plot(T_list, rmse_list, 'o-', ...
    'LineWidth',2, ...
    'Color',[0 0.4470 0.7410]); 
hold on;

plot(T_fit, E_fit, 'r-', ...
    'LineWidth',2);

grid on;

% 坐标轴统一
set(gca, ...
    'FontSize',18, ...
    'LineWidth',1.5, ...
    'TickDir','in', ...
    'Box','on');

xlabel('T (s)', ...
    'FontSize',20, ...
    'FontWeight','bold');

ylabel('RMSE (m)', ...
    'FontSize',20, ...
    'FontWeight','bold');

title('Power-Law Fitting', ...
    'FontSize',21);

legend({'Data','Power-law fit'}, ...
    'FontSize',16, ...
    'Location','best');


% 用拟合模型预测原始点
E_pred = a * T_list.^(-b);

% 计算 R^2
SS_res = sum((rmse_list - E_pred).^2);
SS_tot = sum((rmse_list - mean(rmse_list)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('R^2 = %.6f\n', R2);

RMSE_fit = sqrt(mean((rmse_list - E_pred).^2));
MaxRelErr = max(abs(rmse_list - E_pred)./rmse_list);

fprintf('Prediction RMSE = %.6e\n', RMSE_fit);
fprintf('Max relative error = %.3f%%\n', 100*MaxRelErr);



residual = rmse_list - E_pred;
% 
% figure;
% % subplot(1,2,1)
% plot(T_list, residual, 'o-','LineWidth',1.5);
% grid on;
% xlabel('T (s)');
% ylabel('Residual');
% title('Residual vs T');
% figure;
% % subplot(1,2,2)
% histogram(residual,10);
% grid on;
% title('Residual Histogram');
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
% %% =========图5============
% figure('Color','w');
% clf;
% hold on;
% % ===== 时间筛选 0~1.5 s =====
% idx_time = (t_plot >= 0) & (t_plot <= 1.5);
% t_sel = t_plot(idx_time);
% theta_sel = theta_plot(:, idx_time);
% % θ1 — 实线
% plot(t_sel, theta_sel(1,:), 'b-',  'LineWidth', 2);
% % θ3 — 实线
% plot(t_sel, theta_sel(3,:), 'k-',  'LineWidth', 2);
% % θ2 — 虚线
% plot(t_sel, theta_sel(2,:), 'r--', 'LineWidth', 2);
% grid on;
% xlim([0 1.5]);
% % ===== 坐标轴统一格式 =====
% set(gca, ...
%     'FontSize',18, ...
%     'LineWidth',1.5, ...
%     'TickDir','in', ...
%     'Box','on');
% xlabel('Time (s)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% ylabel('Joint Angle (rad)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% title('Joint Angular Displacements', ...
%     'FontSize',21, ...
%     'FontWeight','normal');
% legend({'$\theta_1$', ...
%         '$\theta_2$', ...
%         '$\theta_3$'}, ...
%         'Interpreter','latex');
% 
% 
% 
% figure('Color','w');
% clf;
% hold on;
% 
% % ===== 时间筛选 =====
% idx_time = (t_plot >= 0) & (t_plot <= 1.5);
% 
% t_sel = t_plot(idx_time);
% dtheta_sel = dtheta_plot(:, idx_time);
% 
% % θ1
% plot(t_sel, dtheta_sel(1,:), 'b-',  'LineWidth', 2);
% % θ3
% plot(t_sel, dtheta_sel(3,:), 'k-',  'LineWidth', 2);
% % θ2
% plot(t_sel, dtheta_sel(2,:), 'r--', 'LineWidth', 2);
% 
% grid on;
% xlim([0 1.5]);
% 
% % 坐标轴风格统一
% set(gca, ...
%     'FontSize',18, ...
%     'LineWidth',1.5, ...
%     'TickDir','in', ...
%     'Box','on');
% 
% xlabel('Time (s)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% 
% ylabel('Angular Velocity (rad/s)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% 
% title('Joint Angular Velocities', ...
%     'FontSize',21);
% 
% legend({'$\dot{\theta}_1$', ...
%         '$\dot{\theta}_2$', ...
%         '$\dot{\theta}_3$'}, ...
%         'Interpreter','latex', ...
%         'FontSize',16);
% 
% figure('Color','w');
% clf;
% hold on;
% 
% % ===== 时间筛选 =====
% idx_time = (t_plot >= 0) & (t_plot <= 1.5);
% 
% t_sel = t_plot(idx_time);
% ddtheta_sel = ddtheta_plot(:, idx_time);
% 
% % θ1
% plot(t_sel, ddtheta_sel(1,:), 'b-',  'LineWidth', 2);
% % θ3
% plot(t_sel, ddtheta_sel(3,:), 'k-',  'LineWidth', 2);
% % θ2
% plot(t_sel, ddtheta_sel(2,:), 'r--', 'LineWidth', 2);
% 
% grid on;
% xlim([0 1.5]);
% 
% % 坐标轴风格统一
% set(gca, ...
%     'FontSize',18, ...
%     'LineWidth',1.5, ...
%     'TickDir','in', ...
%     'Box','on');
% 
% xlabel('Time (s)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% 
% ylabel('Angular Acceleration (rad/s^2)', ...
%     'FontSize',20, ...
%     'FontWeight','bold');
% 
% title('Joint Angular Accelerations', ...
%     'FontSize',21);
% 
% legend({'$\ddot{\theta}_1$', ...
%         '$\ddot{\theta}_2$', ...
%         '$\ddot{\theta}_3$'}, ...
%         'Interpreter','latex', ...
%         'FontSize',16);
% 
