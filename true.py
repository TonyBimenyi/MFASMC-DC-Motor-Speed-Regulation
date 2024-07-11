import numpy as np
import matplotlib.pyplot as plt

# Define parameters
d = 5
rho = 1
eta = 1
lamda = 1
mu = 1
epsilon = 10**-2
alpha = 0.5
T = 0.1
gamma1 = 0.1
gamma2 = 0.1
gamma3 = 0.2
gamma4 = 0.2
rT = 1
m = 1
L = 100

# Generate desired trajectory
yd = np.zeros(L)
for k in range(L):
    yd[k] = 0.5 * np.sin(k * np.pi / 30) + 0.3 * np.cos(k * np.pi / 10)

# Initialize arrays
phi1 = np.zeros(L)
phi2 = np.zeros(L)
phi3 = np.zeros(L)
phi4 = np.zeros(L)

mfa1 = np.zeros(L)
mfa2 = np.zeros(L)
mfa3 = np.zeros(L)
mfa4 = np.zeros(L)

sm1 = np.zeros(L)
sm2 = np.zeros(L)
sm3 = np.zeros(L)
sm4 = np.zeros(L)

u1 = np.zeros(L)
u2 = np.zeros(L)
u3 = np.zeros(L)
u4 = np.zeros(L)

y1 = np.zeros(L + 1)
y2 = np.zeros(L + 1)
y3 = np.zeros(L + 1)
y4 = np.zeros(L + 1)

si1 = np.zeros(L)
si2 = np.zeros(L)
si3 = np.zeros(L)
si4 = np.zeros(L)

# Set initial conditions
y1[0] = 0
y2[0] = 0
y3[0] = 0
y4[0] = 0

# Simulation loop
for k in range(1, L):
    phi1[k] = phi1[k - 1] + eta * u1[k - 1] / (mu + u1[k - 1]**2) * (y1[k] - y1[k - 1] - phi1[k - 1] * u1[k - 1])
    phi2[k] = phi2[k - 1] + eta * u2[k - 1] / (mu + u2[k - 1]**2) * (y2[k] - y2[k - 1] - phi2[k - 1] * u2[k - 1])
    phi3[k] = phi3[k - 1] + eta * u3[k - 1] / (mu + u3[k - 1]**2) * (y3[k] - y3[k - 1] - phi3[k - 1] * u3[k - 1])
    phi4[k] = phi4[k - 1] + eta * u4[k - 1] / (mu + u4[k - 1]**2) * (y4[k] - y4[k - 1] - phi4[k - 1] * u4[k - 1])

    si1[k] = yd[k] - 2 * y1[k] + y4[k]
    si2[k] = y1[k] - 2 * y2[k] + y3[k]
    si3[k] = y2[k] + yd[k] - 2 * y3[k]
    si4[k] = y1[k] + y3[k] - 2 * y4[k]

    mfa1[k] = mfa1[k - 1] + (rho * phi1[k]) / (lamda + abs(phi1[k])**2) * si1[k]
    mfa2[k] = mfa2[k - 1] + (rho * phi2[k]) / (lamda + abs(phi2[k])**2) * si2[k]
    mfa3[k] = mfa3[k - 1] + (rho * phi3[k]) / (lamda + abs(phi3[k])**2) * si3[k]
    mfa4[k] = mfa4[k - 1] + (rho * phi4[k]) / (lamda + abs(phi4[k])**2) * si4[k]

    sm1[k] = mfa1[k] + (yd[k] - y1[k] + alpha * (yd[k] - y1[k]) + epsilon * T * np.sign(k))
    sm2[k] = mfa2[k] + (yd[k] - y2[k] + alpha * (yd[k] - y2[k]) + epsilon * T * np.sign(k))
    sm3[k] = mfa3[k] + (yd[k] - y3[k] + alpha * (yd[k] - y3[k]) + epsilon * T * np.sign(k))
    sm4[k] = mfa4[k] + (yd[k] - y4[k] + alpha * (yd[k] - y4[k]) + epsilon * T * np.sign(k))

    u1[k] = mfa1[k] + gamma1 * sm1[k]
    u2[k] = mfa2[k] + gamma2 * sm2[k]
    u3[k] = mfa3[k] + gamma3 * sm3[k]
    u4[k] = mfa4[k] + gamma4 * sm4[k]

    y1[k + 1] = m / (rT * 8) * u1[k]
    y2[k + 1] = m / (rT * 3) * u2[k]
    y3[k + 1] = m / (rT * 4) * u3[k]
    y4[k + 1] = m / (rT * 5) * u4[k]

# Plot the desired output
plt.figure()
plt.plot(yd, '-b', label='Desired Output')
plt.plot(y1[:-1], '-r', label='Y1')
plt.plot(y2[:-1], '-g', label='Y2')
plt.plot(y3[:-1], '-y', label='Y3')
plt.plot(y4[:-1], '-k', label='Y4')
plt.grid()
plt.legend()
plt.xlabel('Time step')
plt.ylabel('Value')
plt.title('Simulation Results')
plt.show()
