clc
clear
%% 机械臂参数
mm = 1e-3;
L1 = 405*mm;
L_G = 202.5*mm;
L2 = 930*mm;
rb = 218*mm;
rp = 72*mm;

phi = [0;
       2*pi/3;
       4*pi/3];

%% 各点的运动分析
Rz1 = [ cos(phi(1))  -sin(phi(1))  0;
       sin(phi(1))   cos(phi(1))  0
       0            0           1 ];
Rz2 = [ cos(phi(2))  -sin(phi(2))  0;
       sin(phi(2))   cos(phi(2))  0;
       0            0           1 ];

Rz3 = [ cos(phi(3))  -sin(phi(3))  0;
       sin(phi(3))   cos(phi(3))  0;
       0            0           1 ];

syms  x y z vx vy vz ax ay az real

syms theta1 theta2 theta3 real
syms dtheta1 dtheta2 dtheta3 real
syms ddtheta1 ddtheta2 ddtheta3 real

thetai   = [theta1; theta2; theta3];
dthetai  = [dtheta1; dtheta2; dtheta3];
ddthetai = [ddtheta1; ddtheta2; ddtheta3];

v_OG  = sym(zeros(3,3));
v_OB  = sym(zeros(3,3));
a_OG  = sym(zeros(3,3));
a_OB  = sym(zeros(3,3));

for i = 1:3

    Rz = [ cos(phi(i))  -sin(phi(i))  0;
           sin(phi(i))   cos(phi(i))  0;
           0             0            1 ];

    %% -------- 速度 --------
    v_og_local = [ -L_G*dthetai(i)*sin(thetai(i));
                    0;
                   -L_G*dthetai(i)*cos(thetai(i)) ];

    v_ob_local = [ -L1*dthetai(i)*sin(thetai(i));
                    0;
                   -L1*dthetai(i)*cos(thetai(i)) ];

    v_OG(:,i) = Rz * v_og_local;
    v_OB(:,i) = Rz * v_ob_local;

    %% -------- 加速度 --------
    a_og_local = [ -L_G*ddthetai(i)*sin(thetai(i)) - L_G*dthetai(i)^2*cos(thetai(i));
                    0;
                   -L_G*ddthetai(i)*cos(thetai(i)) + L_G*dthetai(i)^2*sin(thetai(i)) ];

    a_ob_local = [ -L1*ddthetai(i)*sin(thetai(i)) - L1*dthetai(i)^2*cos(thetai(i));
                    0;
                   -L1*ddthetai(i)*cos(thetai(i)) + L1*dthetai(i)^2*sin(thetai(i)) ];

    a_OG(:,i) = (Rz * a_og_local);
    a_OB(:,i) = (Rz * a_ob_local);

end
%%  雅各比矩阵推导
v_OCi = [vx; vy; vz];
a_OCi = [ax; ay; az];
a_OP=a_OCi;

I_cell = cell(3,1);
J_cell = cell(3,1);
K_cell = cell(3,1);

for i = 1:3
    
    I_cell{i} = z - L1*sin(thetai(i));
    
    J_cell{i} = x*cos(phi(i)) ...
              + y*sin(phi(i)) ...
              + rp - rb ...
              - L1*cos(thetai(i));
          
    K_cell{i} = y*cos(phi(i)) ...
              - x*sin(phi(i));
end

j_Pi = cell(3,1);

for i = 1:3
    Ii = I_cell{i};
    Ji = J_cell{i};
    Ki = K_cell{i};

    scalar = (Ji*L1*sin(thetai(i)) - Ki*L1*cos(thetai(i)));

    vec = [ 1/(Ji*cos(phi(i)) - Ki*sin(phi(i)));
            1/(Ji*sin(phi(i)) + Ki*cos(phi(i)));
            1/Ii ];

    j_Pi{i} = scalar * vec;   % 3x1
end

syms m1 m2 mp I1 I_motor real
y0 = [0;1;0];
y = cell(3,1);
for i =1:3
        Rz = [ cos(phi(i))  -sin(phi(i))  0;
           sin(phi(i))   cos(phi(i))  0;
           0             0            1 ];
    y{i} = Rz*y0;
end
%% 动力学推导
g = [0;0;-9.8];
tau_p = sym(zeros(3,1));   % 存储 tau_i'
for i = 1:3
    
    Rz = [ cos(phi(i))  -sin(phi(i))  0;
           sin(phi(i))   cos(phi(i))  0;
           0             0            1 ];
    %% --- 末端平台项 ---
    term1 = (mp + 3*m2) * (g - a_OP).' * j_Pi{i};

    %% --- B点项 ---
    vec_B_local = [ -L1*sin(thetai(i));
                     0;
                     L1*cos(thetai(i)) ];

    term2 = m2 * (g - a_OB(:,i)).' * (Rz * vec_B_local);

    %% --- G点项 ---
    vec_G_local = [ -L_G*sin(thetai(i));
                     0;
                     L_G*cos(thetai(i)) ];

    term3 = m1 * (g - a_OG(:,i)).' * (Rz * vec_G_local);

    %% --- 合成 ---
    tau_p(i) = simplify(term1 + term2 + term3);

end

tau = sym(zeros(3,1));
for i = 1:3
tau(i) = (I1+I_motor)*ddthetai(i)-tau_p(i);
end



%% ===== 1. 先提取 G =====
G = subs(tau, [dthetai; ddthetai], zeros(6,1));

%% ===== 2. 剩余惯性+科氏 =====
f1 = tau - G;

%% ===== 3. 提取 M =====
M = jacobian(f1, ddthetai);

%% ===== 4. 剩余科氏 =====
f2 = f1 - M*ddthetai;

%% ===== 5. 提取 C =====
C = sym(zeros(3,3));
for i = 1:3
    for j = 1:3
        C(i,j) = diff(f2(i), dthetai(j));
    end
end
symvar(tau)
tau_func = matlabFunction(tau,'File','tau_delta');
matlabFunction(M,'File','M')
matlabFunction(C,'File','C')
matlabFunction(G,'File','G')































% F = cell(3,1);
% FO = cell(3,1);
% for j = 1:3
% a1=0;
% a2 = 0;
% a3 =0;
% a4 = 0;
% for i = 1:3
%     a1 =a1+ m1*g.'*diff(v_OG(:,i),dthetai(j));
%     a2 = a2 + m2*g.'*diff(v_OB(:,i),dthetai(j));
%     a3 = (mp+3*m2)*g.'*j_Pi{j};
%     a4 = a4 + tau(i)*y{i}.'*y{i}*diff(dthetai(i),dthetai(j));
% end
% F{j}=a1+a2+a3+a4;
% Ftemp = F{j}
% end
% 
% 
% 
% h2 = waitbar(0,'Computing FO...');
% for j = 1:3
% a1=0;
% a2 = 0;
% a3 =0;
% a4 = 0;
% for i = 1:3
%     a1 =a1-m1* a_OG(:,i).'*diff(v_OG(:,i),dthetai(j));
%     a2 = -(mp+3*m2)*a_OP.'*j_Pi{j};
%     a3 = a3-m2*a_OB(:,i).'*diff(v_OB(:,i),dthetai(j));
%     a4 = a4 - (I1+I_motor)*ddthetai(i)*y{i}.'*y{i}*diff(dthetai(i),dthetai(j));
% end
% FO{j}=a1+a2+a3+a4;
% FOtemp = FO{j}
% 
%     waitbar(j/3, h2, sprintf('Computing FO... (%d/3)', j));
% end
% close(h2);
% 
% 
% f= sym(zeros(3,1));
% for i =1:3
%     f(i) = tau(i)-(F{i} + FO{i});
% end
% 
% 
% J = jacobian(f, ddthetai);
% H = jacobian(J(:), ddthetai)

































% % % 提取 M
% M = sym(zeros(3,3));
% for i = 1:3
%     for j = 1:3
%         M(i,j) = diff(f(i), ddthetai(j));
%     end
% end
% check = M.'-M
% % 剩余
% f2 = f - M*ddthetai;
% 
% % 提取 C
% C = sym(zeros(3,3));
% for i = 1:3
%     for j = 1:3
%         C(i,j) = diff(f2(i), dthetai(j));
%     end
% end
% 
% % 提取 G
% G = f2 - C*dthetai;
% 
% % 最终表达式
% tau = M*ddthetai + C*dthetai + G;


















% %% ================= 作图 =================
% figure
% plot(dt, rad2deg(theta(1,:)),'LineWidth',1.2)
% hold on
% plot(dt, rad2deg(theta(2,:)),'LineWidth',1.2)
% plot(dt, rad2deg(theta(3,:)),'LineWidth',1.2)
% grid on
% xlabel('TI_motore (s)')
% ylabel('\theta (deg)')
% title('Joint Angles vs TI_motore')
% legend('\theta_1','\theta_2','\theta_3')
% 
% figure
% plot(dt, rad2deg(dtheta(1,:)),'LineWidth',1.2)
% hold on
% plot(dt, rad2deg(dtheta(2,:)),'LineWidth',1.2)
% plot(dt, rad2deg(dtheta(3,:)),'LineWidth',1.2)
% grid on
% xlabel('TI_motore (s)')
% ylabel('d\theta (deg/s)')
% title('Joint Velocities vs TI_motore')
% legend('d\theta_1','d\theta_2','d\theta_3')
% 
% 
% figure
% plot(dt, rad2deg(ddtheta(1,:)),'LineWidth',1.2)
% hold on
% plot(dt, rad2deg(ddtheta(2,:)),'LineWidth',1.2)
% plot(dt, rad2deg(ddtheta(3,:)),'LineWidth',1.2)
% grid on
% xlabel('TI_motore (s)')
% ylabel('dd\theta (deg/s)')
% title('Joint Acceleration vs TI_motore')
% legend('dd\theta_1','dd\theta_2','dd\theta_3')