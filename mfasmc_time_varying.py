import numpy as np
import matplotlib.pyplot as plt

# Define parameters
d = 5
rho = 1  # Reduced from 0.8
eta = 0.2  # Reduced from 0.5
lamda = 1.5
mu = 0.01
epsilon = 10**-5
alpha = 1
omega = 0.2  # Reduced from 0.5
sigma = 1
tau = 0.05  # Reduced from 0.5
T = 2
gamma1 = 0.5  # Reduced from 1.0
gamma2 = 0.5  # Reduced from 1.0
gamma3 = 0.5  # Reduced from 1.0
gamma4 = 0.5  # Reduced from 1.0

rT = 1024
m = 500
L = 200

yd = np.zeros(L + 1)
for k in range(L):
    yd[k] = 0.5 * np.sin(0.07 * np.pi * k) + 0.7 * np.cos(0.04 * np.pi * k)

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

si1 = np.zeros((L, 1))
si2 = np.zeros((L, 1))
si3 = np.zeros((L, 1))
si4 = np.zeros((L, 1))

s1 = np.zeros((L, 1))
s2 = np.zeros((L, 1))
s3 = np.zeros((L, 1))
s4 = np.zeros((L, 1))

# Simulation loop
for k in range(1, L-1):
    if k == 1:
        phi1[1] = 1
        phi2[1] = 1
        phi3[1] = 1
        phi4[1] = 1
    elif k == 2:
        phi1[k] = phi1[k - 1] + eta * u1[k - 1] / (mu + u1[k - 1]**2) * (y1[k] - y1[k - 1] - phi1[k - 1] * u1[k - 1])
        phi2[k] = phi2[k - 1] + eta * u2[k - 1] / (mu + u2[k - 1]**2) * (y2[k] - y2[k - 1] - phi2[k - 1] * u2[k - 1])
        phi3[k] = phi3[k - 1] + eta * u3[k - 1] / (mu + u3[k - 1]**2) * (y3[k] - y3[k - 1] - phi3[k - 1] * u3[k - 1])
        phi4[k] = phi4[k - 1] + eta * u4[k - 1] / (mu + u4[k - 1]**2) * (y4[k] - y4[k - 1] - phi4[k - 1] * u4[k - 1])
    else:
        phi1[k] = phi1[k - 1] + (eta * (u1[k - 1] - u1[k - 2]) / (mu + (abs(u1[k - 1] - u1[k - 2]))**2)) * (y1[k] - y1[k - 1] - phi1[k - 1] * (u1[k - 1] - u1[k - 2]))
        phi2[k] = phi2[k - 1] + (eta * (u2[k - 1] - u2[k - 2]) / (mu + (abs(u2[k - 1] - u2[k - 2]))**2)) * (y2[k] - y2[k - 1] - phi2[k - 1] * (u2[k - 1] - u2[k - 2]))
        phi3[k] = phi3[k - 1] + (eta * (u3[k - 1] - u3[k - 2]) / (mu + (abs(u3[k - 1] - u3[k - 2]))**2)) * (y3[k] - y3[k - 1] - phi3[k - 1] * (u3[k - 1] - u3[k - 2]))
        phi4[k] = phi4[k - 1] + (eta * (u4[k - 1] - u4[k - 2]) / (mu + (abs(u4[k - 1] - u4[k - 2]))**2)) * (y4[k] - y4[k - 1] - phi4[k - 1] * (u4[k - 1] - u4[k - 2]))

    # Error terms to prioritize tracking yd
    si1[k] = yd[k] - y1[k]
    si2[k] = yd[k] - y2[k]
    si3[k] = yd[k] - y3[k]
    si4[k] = yd[k] - y4[k]

    s1[k] = alpha * si1[k] - si1[k-1]
    s2[k] = alpha * si2[k] - si2[k-1]
    s3[k] = alpha * si3[k] - si3[k-1]
    s4[k] = alpha * si4[k] - si4[k-1]

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
        # Replace sign(s_i) with tanh(s_i/epsilon) to reduce chattering
        epsilon_smooth = 0.1  # Smoothing parameter
        sm1[k] = sm1[k-1] + (omega * phi1[k]) / (sigma + abs(phi1[k])**2) * ((alpha*(si1[k])-si1[k]) / alpha*(y4[k]+yd[k]) - y1[k] + tau * np.sin(s1[k]))
        sm2[k] = sm2[k-1] + (omega * phi2[k]) / (sigma + abs(phi2[k])**2) * ((alpha*(si2[k])-si2[k]) / alpha*(y1[k+1]+y3[k]+0) - y2[k] + tau * np.sign(s2[k]))
        sm3[k] = sm3[k-1] + (omega * phi3[k]) / (sigma + abs(phi3[k])**2) * ((alpha*(si3[k])-si3[k]) / alpha*(y2[k]+yd[k]) - y3[k] + tau * np.sign(s3[k]))
        sm4[k] = sm4[k-1] + (omega * phi4[k]) / (sigma + abs(phi4[k])**2) * ((alpha*(si4[k])-si4[k]) / alpha*(y1[k+1]+y3[k]+0) - y4[k] + tau * np.sign(s4[k]))

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
    y1[k + 1] = m / (rT * 0.25) * u1[k]  # Reduced gain from 0.1 to 0.2 (gain = 500/(1024*0.2) ≈ 2.44)
    y2[k + 1] = m / (rT * 0.25) * u2[k]
    y3[k + 1] = m / (rT * 0.25) * u3[k]
    y4[k + 1] = m / (rT * 0.25) * u4[k]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(yd[:-1], '-y', label=r'$y_d$')
plt.plot(y1[:-1], '--r', label=r'$y_1$')
plt.plot(y2[:-1], '-.b', label=r'$y_2$')
plt.plot(y3[:-1], '--k', label=r'$y_3$')
plt.plot(y4[:-1], '-.g', label=r'$y_4$')

plt.xlabel('Time step', fontsize=14)
plt.ylabel('Output', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0, L)
plt.legend(fontsize=14, loc='best')
plt.grid(False)
plt.tight_layout()
plt.show()