function ddq = dynamics_inverse( ...
    tau_input, q, dq, ...
    I1, I_motor, ...
    ax, ay, az, ...
    m1, m2, mp, ...
    x, y, z)

% 展开
theta1 = q(1);
theta2 = q(2);
theta3 = q(3);

dtheta1 = dq(1);
dtheta2 = dq(2);
dtheta3 = dq(3);

%% ===== 1️⃣ 计算 M =====
M = M_delta( ...
    I1, I_motor, ...
    m1, m2, ...
    theta1, theta2, theta3);

%% ===== 2️⃣ 计算 C =====
C = C_delta( ...
    dtheta1, dtheta2, dtheta3, ...
    m1, m2, ...
    theta1, theta2, theta3);

%% ===== 3️⃣ 计算 G =====
G = G_delta( ...
    ax, ay, az, ...
    m1, m2, mp, ...
    theta1, theta2, theta3, ...
    x, y, z);

%% ===== 4️⃣ 正动力学 =====
ddq = (M + 1e-8*eye(3)) \ (tau_input - C*dq - G);
if any(isnan(M(:)))
    error('M became NaN');
end

if rank(M) < 3
    disp('M rank deficient');
end

condM = cond(M);
if condM > 1e8
    disp(['M badly conditioned, cond = ', num2str(condM)]);
end

end