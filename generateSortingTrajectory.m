function traj = generateSortingTrajectory(A,B,C,D,E,F,G,T,t,nCycle)

if nargin < 10
    nCycle = 1;
end

dt = t(2) - t(1);

%% ===== 生成总时间 =====
t_total = 0:dt:(nCycle*T);
t_mod   = mod(t_total, T);   % 映射到单周期

%% ===== 转列向量 =====
A=A(:); B=B(:); C=C(:); D=D(:);
E=E(:); F=F(:); G=G(:);

% 只取 x,z
A = A([1 3]);
B = B([1 3]);
C = C([1 3]);
D = D([1 3]);
E = E([1 3]);
F = F([1 3]);
G = G([1 3]);

%% =========================
%% ===== 前半段构造 =====
%% =========================

PH_BC = generatePHsegment(B,C);

l_AB = abs(B(2)-A(2));
l_BC = PH_BC.l_all;
l_CD = abs(D(1)-C(1));

l_half = l_AB + l_BC + l_CD;
l_all  = 2*l_half;

%% ===== 时间律 =====
[s,~,~] = timeScaling345(l_all,T,t_mod);

%% ===== 分段 =====
x = zeros(size(t_total));
z = zeros(size(t_total));

S1 = l_AB;
S2 = S1 + l_BC;
S3 = S2 + l_CD;

x_mid = (A(1) + G(1))/2;

for i = 1:length(t_total)

    si = s(i);

    if si <= l_half

        if si <= S1
            ratio = si/l_AB;
            x(i) = A(1);
            z(i) = A(2) + ratio*(B(2)-A(2));

        elseif si <= S2
            gamma = PH_BC.inv_arc(si-S1);
            [x(i),z(i)] = PH_BC.pos(gamma);

        else
            ratio = (si-S2)/l_CD;
            x(i) = C(1) + ratio*(D(1)-C(1));
            z(i) = C(2);
        end

    else

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

        x(i) = 2*x_mid - x_half;
        z(i) = z_half;

    end
end

y = zeros(size(x));

%% ===== 空间速度加速度 =====
vx = gradient(x, dt);
vy = gradient(y, dt);
vz = gradient(z, dt);

ax = gradient(vx, dt);
ay = gradient(vy, dt);
az = gradient(vz, dt);

%% ===== 逆运动学 =====
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

end
