function plotDeltaResults(tVec, tau, traj)
% plotDeltaResults
% 统一绘制 Delta 机器人动力学仿真结果
%
% 输入：
%   tVec  : 时间向量 (1×N)
%   tau   : 3×N 扭矩矩阵
%   traj  : 结构体，包含
%           traj.theta   (3×N)
%           traj.dtheta  (3×N)
%           traj.ddtheta (3×N)

%% ===============================
%  关闭旧窗口
%% ===============================
close all

%% ===============================
%  扭矩
%% ===============================
figure('Name','Joint Torque')
plot(tVec, tau(1,:), 'r','LineWidth',1.5); hold on
plot(tVec, tau(2,:), 'g','LineWidth',1.5);
plot(tVec, tau(3,:), 'b','LineWidth',1.5);

xlabel('Time (s)')
ylabel('\tau (Nm)')
legend('\tau_1','\tau_2','\tau_3')
title('Joint Torque vs Time')
grid on
set(gca,'FontSize',12)

fprintf('---- Torque Info ----\n')
fprintf('Max |tau| = %.2f Nm\n', max(abs(tau),[],'all'))

%% ===============================
%  角度
%% ===============================
figure('Name','Joint Angles')
plot(tVec, rad2deg(traj.theta(1,:)),'LineWidth',1.2); hold on
plot(tVec, rad2deg(traj.theta(2,:)),'LineWidth',1.2)
plot(tVec, rad2deg(traj.theta(3,:)),'LineWidth',1.2)

xlabel('Time (s)')
ylabel('\theta (deg)')
legend('\theta_1','\theta_2','\theta_3')
title('Joint Angles vs Time')
grid on
set(gca,'FontSize',12)

fprintf('Max angle (deg) = %.2f\n', ...
    max(abs(rad2deg(traj.theta)),[],'all'))

%% ===============================
%  角速度
%% ===============================
figure('Name','Joint Velocities')
plot(tVec, rad2deg(traj.dtheta(1,:)),'LineWidth',1.2); hold on
plot(tVec, rad2deg(traj.dtheta(2,:)),'LineWidth',1.2)
plot(tVec, rad2deg(traj.dtheta(3,:)),'LineWidth',1.2)

xlabel('Time (s)')
ylabel('d\theta (deg/s)')
legend('d\theta_1','d\theta_2','d\theta_3')
title('Joint Velocities vs Time')
grid on
set(gca,'FontSize',12)

fprintf('Max velocity (deg/s) = %.2f\n', ...
    max(abs(rad2deg(traj.dtheta)),[],'all'))

%% ===============================
%  角加速度
%% ===============================
figure('Name','Joint Accelerations')
plot(tVec, rad2deg(traj.ddtheta(1,:)),'LineWidth',1.2); hold on
plot(tVec, rad2deg(traj.ddtheta(2,:)),'LineWidth',1.2)
plot(tVec, rad2deg(traj.ddtheta(3,:)),'LineWidth',1.2)

xlabel('Time (s)')
ylabel('dd\theta (deg/s^2)')
legend('dd\theta_1','dd\theta_2','dd\theta_3')
title('Joint Acceleration vs Time')
grid on
set(gca,'FontSize',12)

fprintf('Max acceleration (deg/s^2) = %.2f\n', ...
    max(abs(rad2deg(traj.ddtheta)),[],'all'))

fprintf('----------------------\n\n')

end
