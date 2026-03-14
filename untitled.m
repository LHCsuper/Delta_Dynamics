clc; clear; close all;

%% ===== 控制点 =====
P0 = [0, 0];
P1 = [0, 2];
P3 = [1.5, 3];
P5 = [4, 3];

% 中间控制点（用于构造平滑曲线）
P2 = [0.8, 1.5];
P4 = [2.5, 2.8];

%% ===== 三次Bezier示意曲线 =====
t = linspace(0,1,300)';
B = (1-t).^3 .* P0 + ...
    3*(1-t).^2 .* t .* P2 + ...
    3*(1-t) .* t.^2 .* P4 + ...
    t.^3 .* P5;

%% ===== 绘图 =====
figure; hold on; axis equal;

% 去除刻度和边框
set(gca,'xtick',[],'ytick',[]);
box off;

% 坐标轴
plot([0 4.5],[0 0],'k','LineWidth',1.5);
plot([0 0],[0 3.5],'k','LineWidth',1.5);

% 虚线辅助线
plot([0 P5(1)],[P5(2) P5(2)],'k--','LineWidth',1);
plot([P5(1) P5(1)],[0 P5(2)],'k--','LineWidth',1);

% 曲线
plot(B(:,1), B(:,2),'k','LineWidth',1.8);

% 控制点
plot(P0(1),P0(2),'ko','MarkerFaceColor','k');
plot(P1(1),P1(2),'ko','MarkerFaceColor','k');
plot(P3(1),P3(2),'ko','MarkerFaceColor','k');
plot(P5(1),P5(2),'ko','MarkerFaceColor','k');

% 文字标注（无数字）
% text(P0(1)-0.2,P0(2)-0.2,'P_0(B)');
text(P1(1)-0.6,P1(2),'P_1(P_2)');
text(P3(1)-0.2,P3(2)+0.2,'P_3(P_4)');
text(P5(1)+0.1,P5(2),'P_5(C)');
text(-0.2,3.4,'Y');
text(4.4,-0.2,'X');
text(-0.2,-0.2,'O');

xlim([-0.5 4.5]);
ylim([-0.5 3.5]);
