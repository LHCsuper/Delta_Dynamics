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

T  = 0.4;
dt = 0.001;
t  = 0:dt:T;
nCycle=5;
%% ===== 生成轨迹 =====
traj = generateSortingTrajectory(A,B,C,D,E,F,G,T,t,nCycle);

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

%% PD控制
%% ================= 控制参数 =================
KP = diag([400,400,400]);     % 位置增益
KD = diag([300,300,300]);       % 速度增益

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
        tau(:,k) = KP*e + KD*de;
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

figure;
for i = 1:3
    subplot(3,1,i);

    % 先画 Actual（蓝色实线）
    plot(traj.t, q(i,:), 'b', 'LineWidth', 1.2); 
    hold on;

    % 后画 Desired（红色虚线）→ 在最上层
    plot(traj.t, traj.theta(i,:), 'r--', 'LineWidth', 1.8);

    grid on;
    xlabel('Time (s)');
    ylabel(['\theta_', num2str(i)]);
    legend('Actual','Desired');   % 顺序与绘制一致
end
sgtitle('Delta Joint Tracking Performance');

