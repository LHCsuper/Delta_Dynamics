import numpy as np
import matplotlib.pyplot as plt

# ===== 控制点 =====
P0 = np.array([0.0, 0.0])
P1 = np.array([0.0, 2.0])
P3 = np.array([1.6, 3.0])
P5 = np.array([4.0, 3.0])

# 用Bezier构造平滑示意曲线
P2 = np.array([0.8, 1.4])
P4 = np.array([2.6, 2.8])

t = np.linspace(0, 1, 400)
B = ((1-t)**3)[:, None]*P0 + \
    (3*(1-t)**2*t)[:, None]*P2 + \
    (3*(1-t)*t**2)[:, None]*P4 + \
    (t**3)[:, None]*P5

# ===== 绘图 =====
plt.figure(figsize=(6,5))
plt.gca().set_aspect('equal')

# 去除刻度
plt.xticks([])
plt.yticks([])

# 坐标轴
plt.plot([0, 4.2], [0, 0], 'k', linewidth=1.2)
plt.plot([0, 0], [0, 3.3], 'k', linewidth=1.2)

# 箭头
plt.annotate('', xy=(4.2, 0), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', linewidth=1.2))
plt.annotate('', xy=(0, 3.3), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', linewidth=1.2))

# 虚线
plt.plot([0, P5[0]], [P5[1], P5[1]], 'k--', linewidth=1)
plt.plot([P5[0], P5[0]], [0, P5[1]], 'k--', linewidth=1)

# 曲线
plt.plot(B[:,0], B[:,1], 'k', linewidth=2)

# 控制点
for P in [P0, P1, P3, P5]:
    plt.plot(P[0], P[1], 'ko')

# ===== 标注 =====
plt.text(-0.15, -0.2, r'$O$', fontsize=12)
plt.text(4.25, -0.15, r'$X$', fontsize=12)
plt.text(-0.15, 3.35, r'$Y$', fontsize=12)

plt.text(P0[0]-0.2, P0[1]-0.15, r'$P_0(B)$', fontsize=12)
plt.text(P1[0]-0.7, P1[1], r'$P_1(P_2)$', fontsize=12)
plt.text(P3[0]-0.2, P3[1]+0.2, r'$P_3(P_4)$', fontsize=12)
plt.text(P5[0]+0.1, P5[1], r'$P_5(C)$', fontsize=12)

plt.text(1.7, 1.6, r'$\mathrm{PH\ curve}$', fontsize=13)

# H, r, m 标注
plt.text(-0.35, P5[1]-0.05, r'$H$', fontsize=12)
plt.text((P3[0]+P5[0])/2, P5[1]+0.2, r'$r$', fontsize=12)
plt.text(P5[0]+0.15, P5[1]/2, r'$m$', fontsize=12)

plt.xlim(-0.4, 4.4)
plt.ylim(-0.4, 3.5)
plt.box(False)

plt.tight_layout()
plt.show()
