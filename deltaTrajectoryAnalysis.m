function [v, a, theta, dtheta, ddtheta] = ...
    deltaTrajectoryAnalysis(x, y, z, dt, total_time,plot)

% ===============================
% 输入:
% x,y,z : 轨迹离散点 (Nx1)
% dt : 采样时间
% total_time : 总时间
%
% 输出:
% v      : 笛卡尔速度 (3xN)
% a      : 笛卡尔加速度 (3xN)
% theta  : 关节角 (3xN)
% dtheta : 关节速度 (3xN)
% ddtheta: 关节加速度 (3xN)
% ===============================

N = length(x);

%% ===============================
% 1. 笛卡尔速度
%% ===============================
v = zeros(3,N);
a = zeros(3,N);

pos = [x'; y'; z'];

% ---- 速度 ----
for i = 2:N-1
    v(:,i) = (pos(:,i+1) - pos(:,i-1)) / (2*dt);
end

% 端点
v(:,1)   = (pos(:,2) - pos(:,1)) / dt;
v(:,N)   = (pos(:,N) - pos(:,N-1)) / dt;

% ---- 加速度 ----
for i = 2:N-1
    a(:,i) = (pos(:,i+1) - 2*pos(:,i) + pos(:,i-1)) / (dt^2);
end

a(:,1) = (pos(:,3) - 2*pos(:,2) + pos(:,1)) / (dt^2);
a(:,N) = (pos(:,N) - 2*pos(:,N-1) + pos(:,N-2)) / (dt^2);

%% ===============================
% 2. 逆运动学求 theta
%% ===============================
theta = zeros(3,N);

for i = 1:N
    P = [x(i); y(i); z(i)];
    theta(:,i) = deltaIK(P);
end

%% ===============================
% 3. 关节速度
%% ===============================
dtheta = zeros(3,N);
ddtheta = zeros(3,N);

% ---- 关节速度 ----
for i = 2:N-1
    dtheta(:,i) = (theta(:,i+1) - theta(:,i-1)) / (2*dt);
end

dtheta(:,1) = (theta(:,2) - theta(:,1)) / dt;
dtheta(:,N) = (theta(:,N) - theta(:,N-1)) / dt;

% ---- 关节加速度 ----
for i = 2:N-1
    ddtheta(:,i) = (theta(:,i+1) - 2*theta(:,i) + theta(:,i-1)) / (dt^2);
end

ddtheta(:,1) = (theta(:,3) - 2*theta(:,2) + theta(:,1)) / (dt^2);
ddtheta(:,N) = (theta(:,N) - 2*theta(:,N-1) + theta(:,N-2)) / (dt^2);


end


