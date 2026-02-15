function traj = generateDeltaTrajectory(trajFunc, tVec, deltaIK)

n = length(tVec);
dt = tVec(2) - tVec(1);

%% 预分配
theta  = zeros(3,n);
dtheta = zeros(3,n);
ddtheta = zeros(3,n);

x = zeros(1,n);
y = zeros(1,n);
z = zeros(1,n);

vx = zeros(1,n);
vy = zeros(1,n);
vz = zeros(1,n);

ax = zeros(1,n);
ay = zeros(1,n);
az = zeros(1,n);

%% 生成末端轨迹 + IK
for k = 1:n
    
    t = tVec(k);
    
    P = trajFunc(t);
    
    x(k) = P(1);
    y(k) = P(2);
    z(k) = P(3);
    
    theta(:,k) = deltaIK(P);  % 必须是弧度
    
end

%% 数值微分（中心差分）
for i = 1:3
    
    dtheta(i,2:end-1) = (theta(i,3:end) - theta(i,1:end-2)) / (2*dt);
    
    ddtheta(i,2:end-1) = (theta(i,3:end) ...
                        - 2*theta(i,2:end-1) ...
                        + theta(i,1:end-2)) / (dt^2);
end

%% 末端速度加速度
vx(2:end-1) = (x(3:end) - x(1:end-2)) / (2*dt);
vy(2:end-1) = (y(3:end) - y(1:end-2)) / (2*dt);
vz(2:end-1) = (z(3:end) - z(1:end-2)) / (2*dt);

ax(2:end-1) = (x(3:end) - 2*x(2:end-1) + x(1:end-2)) / (dt^2);
ay(2:end-1) = (y(3:end) - 2*y(2:end-1) + y(1:end-2)) / (dt^2);
az(2:end-1) = (z(3:end) - 2*z(2:end-1) + z(1:end-2)) / (dt^2);

%% =========================
%% ======= 调试打印 =======
%% =========================

fprintf('\n========= Trajectory Debug Info =========\n');

maxTheta = max(abs(theta(:)));
maxDtheta = max(abs(dtheta(:)));
maxDDtheta = max(abs(ddtheta(:)));

maxPos = max(abs([x y z]));
maxAcc = max(abs([ax ay az]));

fprintf('Max |theta|        = %.6f rad (%.2f deg)\n', ...
        maxTheta, rad2deg(maxTheta));

fprintf('Max |dtheta|       = %.6f rad/s\n', maxDtheta);
fprintf('Max |ddtheta|      = %.6f rad/s^2\n', maxDDtheta);

fprintf('Max |position|     = %.6f m\n', maxPos);
fprintf('Max |acceleration| = %.6f m/s^2\n', maxAcc);

% 自动单位判断
if maxTheta > 2*pi
    warning('theta 可能是度而不是弧度！');
end

if maxPos > 5
    warning('位置单位可能是 mm 而不是 m！');
end


%% 打包输出
traj.theta  = theta;
traj.dtheta = dtheta;
traj.ddtheta = ddtheta;

traj.x = x;
traj.y = y;
traj.z = z;

traj.vx = vx;
traj.vy = vy;
traj.vz = vz;

traj.ax = ax;
traj.ay = ay;
traj.az = az;

traj.t = tVec;
%% =========================
%% ======= 绘图输出 =======
%% =========================

figure

subplot(3,3,1)
plot(traj.t, traj.x,'LineWidth',1.2)
title('x(t)')
xlabel('Time (s)')
ylabel('x (m)')
grid on

subplot(3,3,2)
plot(traj.t, traj.y,'LineWidth',1.2)
title('y(t)')
xlabel('Time (s)')
ylabel('y (m)')
grid on

subplot(3,3,3)
plot(traj.t, traj.z,'LineWidth',1.2)
title('z(t)')
xlabel('Time (s)')
ylabel('z (m)')
grid on

subplot(3,3,4)
plot(traj.t, traj.vx,'LineWidth',1.2)
title('vx(t)')
xlabel('Time (s)')
ylabel('vx (m/s)')
grid on

subplot(3,3,5)
plot(traj.t, traj.vy,'LineWidth',1.2)
title('vy(t)')
xlabel('Time (s)')
ylabel('vy (m/s)')
grid on

subplot(3,3,6)
plot(traj.t, traj.vz,'LineWidth',1.2)
title('vz(t)')
xlabel('Time (s)')
ylabel('vz (m/s)')
grid on

subplot(3,3,7)
plot(traj.t, traj.ax,'LineWidth',1.2)
title('ax(t)')
xlabel('Time (s)')
ylabel('ax (m/s^2)')
grid on

subplot(3,3,8)
plot(traj.t, traj.ay,'LineWidth',1.2)
title('ay(t)')
xlabel('Time (s)')
ylabel('ay (m/s^2)')
grid on

subplot(3,3,9)
plot(traj.t, traj.az,'LineWidth',1.2)
title('az(t)')
xlabel('Time (s)')
ylabel('az (m/s^2)')
grid on
end
