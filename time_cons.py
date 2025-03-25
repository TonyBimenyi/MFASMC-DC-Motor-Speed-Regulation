import numpy as np
import matplotlib.pyplot as plt

# Define parameters
d = 5
rho = 1
eta = 1
lamda = 50
mu = 1
epsilon = 10**-5
alpha = 1
T = 0.1
gamma1 = 0.15
gamma2 = 0.15
gamma3 = 0.45
gamma4 = 0.45
rT = 1024
m = 350
L = 200
font_style = 'Times New Roman'

# Generate constant desired trajectory
yd = 0.5 * (1 + np.sign(np.sin(np.linspace(0, 2 * np.pi, L + 1))))  # Square wave

# yd = np.zeros(L + 1)
# for k in range(L):
#     yd[k] = 0.5 * np.sin(0.07 * np.pi * k) + 0.7 * np.cos(0.04 * np.pi * k)



# Initialize arrays
phi1 = np.zeros((L, 1))
phi2 = np.zeros((L, 1))
phi3 = np.zeros((L, 1))
phi4 = np.zeros((L, 1))

mfa1 = np.zeros((L, 1))
mfa2 = np.zeros((L, 1))
mfa3 = np.zeros((L, 1))
mfa4 = np.zeros((L, 1))

sm1 = np.zeros((L, 1))
sm2 = np.zeros((L, 1))
sm3 = np.zeros((L, 1))
sm4 = np.zeros((L, 1))

u1 = np.zeros((L, 1))
u2 = np.zeros((L, 1))
u3 = np.zeros((L, 1))
u4 = np.zeros((L, 1))

y1 = np.zeros((L + 1, 1))
y2 = np.zeros((L + 1, 1))
y3 = np.zeros((L + 1, 1))
y4 = np.zeros((L + 1, 1))

e1 = np.zeros((L + 1, 1))
e2 = np.zeros((L + 1, 1))
e3 = np.zeros((L + 1, 1))
e4 = np.zeros((L + 1, 1))

si1 = np.zeros((L, 1))
si2 = np.zeros((L, 1))
si3 = np.zeros((L, 1))
si4 = np.zeros((L, 1))



# Simulation loop
for k in range(1, L-1):
    if k == 0:
        phi1[0] = 1
        phi2[0] = 1
        phi3[0] = 1
        phi4[0] = 1
    elif k == 1:
        phi1[k] = phi1[k - 1] + eta * u1[k - 1] / (mu + u1[k - 1]**2) * (y1[k] - y1[k - 1] - phi1[k - 1] * u1[k - 1])
        phi2[k] = phi2[k - 1] + eta * u2[k - 1] / (mu + u2[k - 1]**2) * (y2[k] - y2[k - 1] - phi2[k - 1] * u2[k - 1])
        phi3[k] = phi3[k - 1] + eta * u3[k - 1] / (mu + u3[k - 1]**2) * (y3[k] - y3[k - 1] - phi3[k - 1] * u3[k - 1])
        phi4[k] = phi4[k - 1] + eta * u4[k - 1] / (mu + u4[k - 1]**2) * (y4[k] - y4[k - 1] - phi4[k - 1] * u4[k - 1])
    else:
        phi1[k] = phi1[k - 1] + (eta * (u1[k - 1] - u1[k - 2]) / (mu + (abs(u1[k - 1] - u1[k - 2]))**2)) * (y1[k] - y1[k - 1] - phi1[k - 1] * (u1[k - 1] - u1[k - 2]))
        phi2[k] = phi2[k - 1] + (eta * (u2[k - 1] - u2[k - 2]) / (mu + (abs(u2[k - 1] - u2[k - 2]))**2)) * (y2[k] - y2[k - 1] - phi2[k - 1] * (u2[k - 1] - u2[k - 2]))
        phi3[k] = phi3[k - 1] + (eta * (u3[k - 1] - u3[k - 2]) / (mu + (abs(u3[k - 1] - u3[k - 2]))**2)) * (y3[k] - y3[k - 1] - phi3[k - 1] * (u3[k - 1] - u3[k - 2]))
        phi4[k] = phi4[k - 1] + (eta * (u4[k - 1] - u4[k - 2]) / (mu + (abs(u4[k - 1] - u4[k - 2]))**2)) * (y4[k] - y4[k - 1] - phi4[k - 1] * (u4[k - 1] - u4[k - 2]))

    si1[k] = yd[k] - 2 * y1[k] + y4[k]
    si2[k] = y1[k] - 2 * y2[k] + y3[k]
    si3[k] = y2[k] + yd[k] - 2 * y3[k]
    si4[k] = y1[k] + y3[k] - 2 * y4[k]

    if k == 1:
        mfa1[0] = 0
        mfa2[0] = 0
        mfa3[0] = 0
        mfa4[0] = 0
    else:
        mfa1[k] = mfa1[k - 1] + (rho * phi1[k]) / (lamda + abs(phi1[k])**2) * si1[k]
        mfa2[k] = mfa2[k - 1] + (rho * phi2[k]) / (lamda + abs(phi2[k])**2) * si2[k]
        mfa3[k] = mfa3[k - 1] + (rho * phi3[k]) / (lamda + abs(phi3[k])**2) * si3[k]
        mfa4[k] = mfa4[k - 1] + (rho * phi4[k]) / (lamda + abs(phi4[k])**2) * si4[k]

    if k == 1:
        sm1[0] = 0
        sm2[0] = 0
        sm3[0] = 0
        sm4[0] = 0
    else:
        sm1[k] = sm1[k - 1] + (yd[k + 1] - y1[k] + alpha * yd[k + 1] - y1[k] + epsilon * T * np.sign(k))
        sm2[k] = sm2[k - 1] + (yd[k + 1] - y1[k] + alpha * yd[k + 1] - y1[k] + epsilon * T * np.sign(k))
        sm3[k] = sm3[k - 1] + (yd[k + 1] - y1[k] + alpha * yd[k + 1] - y1[k] + epsilon * T * np.sign(k))
        sm4[k] = sm4[k - 1] + (yd[k + 1] - y1[k] + alpha * yd[k + 1] - y1[k] + epsilon * T * np.sign(k))

    if k == 1:
        u1[0] = 0.1
        u2[0] = 0.1
        u3[0] = 0.1
        u4[0] = 0.1
    else:
        u1[k] = mfa1[k] + gamma1 * sm1[k]
        u2[k] = mfa2[k] + gamma2 * sm2[k]
        u3[k] = mfa3[k] + gamma3 * sm3[k]
        u4[k] = mfa4[k] + gamma4 * sm4[k]

    y1[0] = 0.1
    y2[0] = 0.1
    y3[0] = 0.1
    y4[0] = 0.1
    y1[k + 1] = m / (rT * 0.1) * u1[k]
    y2[k + 1] = m / (rT * 0.1) * u2[k]
    y3[k + 1] = m / (rT * 0.3) * u3[k]
    y4[k + 1] = m / (rT * 0.3) * u4[k]

    e1[k] = yd[k] - y1[k]
    e2[k] = yd[k] - y2[k]
    e3[k] = yd[k] - y3[k]
    e4[k] = yd[k] - y4[k]

# Create subplots
# fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# # Plot each output and its tracking error
# axs[0, 0].plot(y1[:-1], '-*r', markersize=7, label='y1')
# axs[0, 0].plot(yd, '-b', label='yd')
# # axs[0, 0].set_title('Tracking Performance y1')
# axs[0, 0].set_xlim(0, L)
# axs[0, 0].set_ylim(min(y1.min(), yd.min()) - 0.1, max(y1.max(), yd.max()) + 0.1)
# axs[0, 0].legend()
# axs[0, 0].grid()

# axs[0, 1].plot(y2[:-1], '-xg', markersize=7, label='y2')
# axs[0, 1].plot(yd, '-b', label='yd')
# # axs[0, 1].set_title('Tracking Performance y2')
# axs[0, 1].set_xlim(0, L)
# axs[0, 1].set_ylim(min(y2.min(), yd.min()) - 0.1, max(y2.max(), yd.max()) + 0.1)
# axs[0, 1].legend()
# axs[0, 1].grid()

# axs[1, 0].plot(y3[:-1], '-oy', label='y3', markersize=5)
# axs[1, 0].plot(yd, '-b', label='yd')
# # axs[1, 0].set_title('Tracking Performance y3')
# axs[1, 0].set_xlim(0, L)
# axs[1, 0].set_ylim(min(y3.min(), yd.min()) - 0.1, max(y3.max(), yd.max()) + 0.1)
# axs[1, 0].legend()
# axs[1, 0].grid()

# axs[1, 1].plot(y4[:-1], '-vk', label='y4', markersize=7)
# axs[1, 1].plot(yd, '-b', label='yd')
# # axs[1, 1].set_title('Tracking Performance y4')
# axs[1, 1].set_xlim(0, L)
# # axs[1, 1].set_ylim(-1.3, 1.3)
# axs[1, 1].set_ylim(min(y4.min(), yd.min()) - 0.1, max(y4.max(), yd.max()) + 0.1)
# axs[1, 1].legend()
# axs[1, 1].grid()

# # Set common labels
# for ax in axs.flat:
#     ax.set_xlabel('Step')
#     ax.set_ylabel('Output')

# plt.tight_layout()
# plt.show()

plt.figure(figsize=(10, 6))  # Adjust figure size for better visibility
plt.rcParams['font.family'] = 'Times New Roman'  # Set font globally to Times New Roman

plt.plot(yd[:-1], '-y', linewidth=2.5, label=r'$y_d$')
plt.plot(y1[:-1], '--r', linewidth=2.5, label=r'$y_1$')
plt.plot(y2[:-1], '-.b', linewidth=2.5, label=r'$y_2$')
plt.plot(y3[:-1], '--k', linewidth=2.5, label=r'$y_3$')
plt.plot(y4[:-1], '-.g', linewidth=2.5, label=r'$y_4$')

# Enlarge font size for axis labels and legend
plt.xlabel('Time Step', fontsize=16, fontname=font_style)
plt.ylabel('Output', fontsize=16, fontname=font_style)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlim(0,L)  # Set y-axis to range from -0.35 to 0.35

# Improve legend readability
plt.legend(fontsize=16, loc='best')
# plt.title('Tracking performance',fontsize=15, fontweight='bold',fontname=font_style)

plt.grid(False)  # Add grid for better readability
plt.tight_layout()  # Adjust layout for clarity

plt.show()