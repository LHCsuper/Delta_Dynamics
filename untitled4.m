clc; clear; close all;

%% ===== 定义关键点（X,Z）=====
A = [-0.3, -0.7];   % 左下
B = [-0.3, -0.6];   % 左上
C = [ 0.3, -0.6];   % 右上
D = [ 0.3, -0.7];   % 右下

%% ===== 轨迹拼接 =====
X = [A(1) B(1) C(1) D(1)];
Z = [A(2) B(2) C(2) D(2)];

%% 图2
figure('Color','w');   % 白底

plot(X, Z, ...
    'LineWidth', 2, ...
    'Color',[0 0.4470 0.7410]);   % 默认MATLAB蓝
hold on;

% 关键点
scatter(X, Z, 80, ...
    'filled', ...
    'MarkerFaceColor',[0.85 0.33 0.1], ...
    'MarkerEdgeColor','k');

grid on;
axis equal;

% 坐标范围
xlim([-0.35 0.35]);
ylim([-0.9 -0.4]);

% 坐标轴统一格式
set(gca, ...
    'FontSize',18, ...
    'LineWidth',1.5, ...
    'TickDir','in', ...
    'Box','on');

xlabel('X (m)', ...
    'FontSize',20, ...
    'FontWeight','bold');

ylabel('Z (m)', ...
    'FontSize',20, ...
    'FontWeight','bold');

title('Simple Gate Shape in X-Z Plane', ...
    'FontSize',21, ...
    'FontWeight','normal');