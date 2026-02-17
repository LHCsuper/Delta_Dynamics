function PH = generatePHsegment(P0, P5)

P0 = P0(:);
P5 = P5(:);

dx = P5(1) - P0(1);
dz = P5(2) - P0(2);

r = abs(dx);
m = abs(dz);

sx = sign(dx); if sx==0, sx=1; end
sz = sign(dz); if sz==0, sz=1; end

alpha = sqrt(m^2 + r^2 + 70*m*r);

u0 = sqrt(5*(35*m + r - alpha)/34);
u2 = sqrt(5*(m + 35*r - alpha)/68);
v2 = u2;

% 总弧长
l_all = (1/5)*u0^2 + (1/15)*u0*u2 + (1/5)*u2^2 + (1/5)*v2^2;

% ===== 弧长查表 =====
N = 4000;
gamma_tab = linspace(0,1,N);

coefA = (u0^2 + 2*u0*u2 + u2^2 + v2^2);
coefB = (4*u0^2 + 4*u0*u2);
coefC = (6*u0^2 + 2*u0*u2);

l_tab = coefA*(gamma_tab.^5)/5 ...
      - coefB*(gamma_tab.^4)/4 ...
      + coefC*(gamma_tab.^3)/3 ...
      - 4*u0^2*(gamma_tab.^2)/2 ...
      + u0^2*gamma_tab;

% ===== 输出结构 =====
PH.l_all = l_all;

PH.inv_arc = @(s) interp1(l_tab, gamma_tab, s, 'pchip');

PH.pos = @(g) localPos(g, P0, u0, u2, v2, sx, sz);

end


% ===== 局部函数 =====
function [x,z] = localPos(g, P0, u0, u2, v2, sx, sz)

x_local = (2/5)*(u0*v2 + u2*v2).*g.^5 ...
        - (u0*v2).*g.^4 ...
        + (2/3)*(u0*v2).*g.^3;

y_local = u0^2.*g.*(1-g).^4 ...
        + 2*u0^2.*g.^2.*(1-g).^3 ...
        + (2*u0^2 + (2/3)*u0*u2).*g.^3.*(1-g).^2 ...
        + (u0^2 + (1/3)*u0*u2).*g.^4.*(1-g) ...
        + ((1/5)*u0^2 + (1/15)*u0*u2).*g.^5;

x = P0(1) + sx*x_local;
z = P0(2) + sz*y_local;

end
