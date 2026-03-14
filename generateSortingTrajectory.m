function traj = generateSortingTrajectory(A,B,C,D,E,F,G,T,t,nCycle)
% 一个闭合周期 = A -> G -> A
% 其中 T 表示单程 A->G 的时间；闭合周期时间为 2T
% nCycle 表示闭合周期个数

if nargin < 10
    nCycle = 1;
end

dt = t(2) - t(1);

%% ===== 转列向量 & 只取 x,z =====
A=A(:); B=B(:); C=C(:); D=D(:);
E=E(:); F=F(:); G=G(:);

A = A([1 3]);
B = B([1 3]);
C = C([1 3]);
D = D([1 3]);
E = E([1 3]);
F = F([1 3]);
G = G([1 3]);

%% =========================
%% ===== 先生成“单程 A->G”轨迹（时长 T）=====
%% =========================

% ---- PH 段（BC） ----
PH_BC = generatePHsegment(B,C);

l_AB = abs(B(2)-A(2));
l_BC = PH_BC.l_all;
l_CD = abs(D(1)-C(1));

l_half = l_AB + l_BC + l_CD;   % A->D
l_all  = 2*l_half;             % A->G（镜像得到）

% ---- 单程时间轴：不包含端点 T，避免与回程拼接时重复G ----
t_fwd = 0:dt:(T - dt);         % [0, T)
t_mod = t_fwd;                 % 单程不需要 mod（也更稳健）

% ---- 单程时间律：s 从 0 -> l_all ----
[s_fwd,~,~] = timeScaling345(l_all, T, t_mod);

% ---- 单程分段求 x,z（A->G）----
x_fwd = zeros(size(t_fwd));
z_fwd = zeros(size(t_fwd));

S1 = l_AB;
S2 = S1 + l_BC;
S3 = S2 + l_CD; %#ok<NASGU>  % 预留，当前用不到 S3

x_mid = (A(1) + G(1))/2;

for i = 1:length(t_fwd)

    si = s_fwd(i);

    if si <= l_half
        % ===== A->D =====
        if si <= S1
            ratio = si/l_AB;
            x_fwd(i) = A(1);
            z_fwd(i) = A(2) + ratio*(B(2)-A(2));

        elseif si <= S2
            gamma = PH_BC.inv_arc(si-S1);
            [x_fwd(i),z_fwd(i)] = PH_BC.pos(gamma);

        else
            ratio = (si-S2)/l_CD;
            x_fwd(i) = C(1) + ratio*(D(1)-C(1));
            z_fwd(i) = C(2);
        end

    else
        % ===== D->G（镜像）=====
        s_mirror = l_all - si;

        if s_mirror <= S1
            ratio = s_mirror/l_AB;
            x_half = A(1);
            z_half = A(2) + ratio*(B(2)-A(2));

        elseif s_mirror <= S2
            gamma = PH_BC.inv_arc(s_mirror-S1);
            [x_half,z_half] = PH_BC.pos(gamma);

        else
            ratio = (s_mirror-S2)/l_CD;
            x_half = C(1) + ratio*(D(1)-C(1));
            z_half = C(2);
        end

        x_fwd(i) = 2*x_mid - x_half;
        z_fwd(i) = z_half;
    end
end

% 单程 y=0
y_fwd = zeros(size(x_fwd));

%% =========================
%% ===== 构造“闭合周期 A->G->A”（时长 2T）=====
%% =========================

% 回程：G->A（反向单程），去掉首点避免重复 G
x_bwd = fliplr(x_fwd);
y_bwd = fliplr(y_fwd);
z_bwd = fliplr(z_fwd);

% 拼接得到一个闭合周期（长度 N_cycle）
x_cycle = [x_fwd, x_bwd];
y_cycle = [y_fwd, y_bwd];
z_cycle = [z_fwd, z_bwd];

% s 也拼一下（用于调试/可视化）：回程可用 l_all -> 0
s_bwd   = fliplr(s_fwd);
s_cycle = [s_fwd, s_bwd];

% 闭合周期时间轴（不包含 2T 端点）
T_cycle = 2*T;
t_cycle = 0:dt:(T_cycle - dt);

% 保证长度一致（保险）
N_cycle = length(t_cycle);
x_cycle = x_cycle(1:N_cycle);
y_cycle = y_cycle(1:N_cycle);
z_cycle = z_cycle(1:N_cycle);
s_cycle = s_cycle(1:N_cycle);

%% =========================
%% ===== 扩展到 nCycle 个闭合周期 =====
%% =========================

if nCycle == 1
    x = x_cycle; y = y_cycle; z = z_cycle;
    s = s_cycle;
    t_total = t_cycle;
else
    % 重复拼接：每次去掉首点避免重复周期边界点
    x = x_cycle; y = y_cycle; z = z_cycle; s = s_cycle;
    for k = 2:nCycle
        x = [x, x_cycle(2:end)];
        y = [y, y_cycle(2:end)];
        z = [z, z_cycle(2:end)];
        s = [s, s_cycle(2:end)];
    end
    t_total = 0:dt:(length(x)-1)*dt;
end

%% ===== 空间速度加速度（基于最终 x,y,z）=====
vx = gradient(x, dt);
vy = gradient(y, dt);
vz = gradient(z, dt);

ax = gradient(vx, dt);
ay = gradient(vy, dt);
az = gradient(vz, dt);

%% ===== 逆运动学（基于最终轨迹长度）=====
N = length(t_total);
theta = zeros(3,N);

for i = 1:N
    P = [x(i); y(i); z(i)];
    theta(:,i) = deltaIK(P);
end

%% ===== 关节速度加速度 =====
dtheta  = gradient(theta, dt);
ddtheta = gradient(dtheta, dt);

%% ===== 输出 =====
traj.x = x;
traj.y = y;
traj.z = z;

traj.v = [vx; vy; vz];
traj.a = [ax; ay; az];

traj.theta   = theta;
traj.dtheta  = dtheta;
traj.ddtheta = ddtheta;

traj.s = s;
traj.t = t_total;

traj.T_single = T;        % 单程时间
traj.T_cycle  = T_cycle;  % 闭合周期时间

end
