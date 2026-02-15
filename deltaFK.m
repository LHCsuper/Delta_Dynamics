function P = deltaFK(theta)
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

J11 = rp - rb - L1*cos(theta(1));
J12 = rp - rb - L1*cos(theta(2));
J13 = rp - rb - L1*cos(theta(3));

I11 = J11^2 + (L1*sin(theta(1)))^2;
I12 = J12^2 + (L1*sin(theta(2)))^2;
I13 = J13^2 + (L1*sin(theta(3)))^2;
c1 = cos(phi(1));  s1 = sin(phi(1));
c2 = cos(phi(2));  s2 = sin(phi(2));
c3 = cos(phi(3));  s3 = sin(phi(3));

% ---------- 差分项 ----------
% ===== 由(7-25)(7-26)直接构造线性方程：M*[x;y] = -g*z - d =====

alpha12 = 2*(J11*c1 - J12*c2);
beta12  = 2*(J11*s1 - J12*s2);
gamma12 = 2*L1*(sin(theta(1)) - sin(theta(2)));
delta12 = (I11 - I12);

alpha13 = 2*(J11*c1 - J13*c3);
beta13  = 2*(J11*s1 - J13*s3);
gamma13 = 2*L1*(sin(theta(1)) - sin(theta(3)));
delta13 = (I11 - I13);

M = [alpha12, beta12;
     alpha13, beta13];

tol = 1e-12;
if abs(det(M)) < tol
    error('Delta FK singular configuration: det(M) ~ 0');
end

g = [gamma12; gamma13];
d = [delta12; delta13];

% x,y = x0 + x1*z, y0 + y1*z
xy0 = -(M \ d);   % 常数项
xy1 = -(M \ g);   % z 系数

D1 = xy0(1);
D3 = xy0(2);
D2 = xy1(1);
D4 = xy1(2);



% ===== 直接选用 i = 1 =====
i = 1;

ci = cos(phi(i));
si = sin(phi(i));

J1i = [J11,J12,J13];
I1i = [I11,I12,I13];

a = D2^2 + D4^2 + 1;
b = D1*D2 + D3*D4 + J1i(i)*(D2*ci + D4*si) + L1*sin(theta(i));
c = D1^2 + D3^2 + 2*J1i(i)*(D1*ci + D3*si) + I1i(i) - L2^2;

Delta = b^2 - a*c;

if Delta < 0
    error('Delta FK: no real solution');
end

z = (-b - sqrt(Delta))/a;   % 选向下根
x = D1 + D2*z;
y = D3 + D4*z;
P = [x; y; z];


end