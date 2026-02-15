function theta = deltaIK(P)

% ===============================
% Delta Robot Inverse Kinematics
% Input:
%   P = [x; y; z]   (单位 mm)
%
% Output:
%   theta = [θ1; θ2; θ3]  (单位 rad)
% ===============================
mm = 1e-3;
% -------- 参数 --------
L1 = 405*mm;
L2 = 930*mm;
rb = 218*mm;
rp = 72*mm;

phi = [0;
       2*pi/3;
       4*pi/3];

% -------- 末端位置 --------
x = P(1);
y = P(2);
z = P(3);

isSymMode = isa(P,'sym');

if isSymMode
    theta = sym(zeros(3,1));
else
    theta = zeros(3,1);
end

for i = 1:3
    
    % 局部旋转坐标
    Xi = x*cos(phi(i)) + y*sin(phi(i));
    Yi = -x*sin(phi(i)) + y*cos(phi(i));
    
    Di = Xi + rp - rb;
    
    % Ki
    Ki = (Di^2 + Yi^2 + z^2 + L1^2 - L2^2)/(2*L1);
    
    % Ri
    Ri = sqrt(Di^2 + z^2);
    
    % 可达性检查
    if ~isa(Ki,'sym') && ~isa(Ri,'sym')
        if abs(Ki/Ri) > 1
            error('该位置超出工作空间');
        end
    end
    
    % alpha
    alpha = atan2(z, Di);
    
    % 选肘下解（负号）
    theta(i) = alpha + acos(Ki/Ri);
    
end
if isSymMode
    theta = simplify(theta);
end
theta = -theta;
end

