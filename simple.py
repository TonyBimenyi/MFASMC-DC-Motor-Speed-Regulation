import numpy as np
import matplotlib.pyplot as plt

# Define parameters
d = 5
rho = 1
eta = 1
lamda = 20
mu = 1
epsilon = 10**-5
alpha = 1  # Adjusted
T = 0.1
gamma1 = 0.15  # Increased
gamma2 = 0.15 # Increased
gamma3 = 0.45
gamma4 = 0.45
rT = 1024
m = 350
L = 200


# Generate desired trajectory
# yd = np.zeros(L)
# for k in range(L):
#     yd[k] = 0.4 * np.sin(0.1 * k) + 0.3 * np.sin(0.05 * k) + 0.2 * np.sin(0.01 * k)

    # Generate a more complex desired trajectory
# yd = np.zeros(L + 1)
# for k in range(L):
#     yd[k] = 0.5 * np.sin(k * np.pi / 30) + 0.3 * np.cos(k * np.pi / 10)


# Generate constant desired trajectory
yd = 0.5 * (1 + np.sign(np.sin(np.linspace(0, 2 * np.pi, L + 1))))  # Square wave

# # Generate desired trajectory
# yd = np.zeros(L + 1)
# for k in range(L):
#     yd[k] = 0.5 * np.sin(k * np.pi / 30) + 0.3 * np.cos(k * np.pi / 10)
    

# Initialize arrays
phi1 = np.zeros((L,1))
phi2 = np.zeros((L,1))
phi3 = np.zeros((L,1))
phi4 = np.zeros((L,1))

mfa1 = np.zeros((L,1))
mfa2 = np.zeros((L,1))
mfa3 = np.zeros((L,1))
mfa4 = np.zeros((L,1))

sm1 = np.zeros((L,1))
sm2 = np.zeros((L,1))
sm3 = np.zeros((L,1))
sm4 = np.zeros((L,1))

u1 = np.zeros((L,1))
u2 = np.zeros((L,1))
u3 = np.zeros((L,1))
u4 = np.zeros((L,1))

y1 = np.zeros((L + 1,1))
y2 = np.zeros((L + 1,1))
y3 = np.zeros((L + 1,1))
y4 = np.zeros((L + 1,1))

e1 = np.zeros((L + 1,1))
e2 = np.zeros((L + 1,1))
e3 = np.zeros((L + 1,1))
e4 = np.zeros((L + 1,1))

si1 = np.zeros((L,1))
si2 = np.zeros((L,1))
si3 = np.zeros((L,1))
si4 = np.zeros((L,1))




# Simulation loop
for k in range(1,L-1):
    if k == 0:
        phi1[0] = 1  # the initial value of phi(k), can't set to zero.
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

# Plot the desired output
plt.figure()
# plt.plot(yd, '-b', label='Desired Output')
# plt.plot(y1[:-1], '-*r', markersize=4, label='Y1')
# plt.plot(y2[:-1], '-og', markersize=4, label='Y2')
# plt.plot(y3[:-1], '--y', label='Y3')
# plt.plot(y4[:-1], '-k', label='Y4')

plt.plot(phi1[:-1], '-*r', markersize=4, label='phi1')
plt.plot(phi2[:-1], '-+g', markersize=4, label='phi2')
plt.plot(phi3[:-1], '--y', label='phi3')
plt.plot(phi4[:-1], '-b', label='phi4')
plt.grid()
plt.legend()


# Set axis limits
plt.xlim(0, L)  # Set x-axis to range from 0 to 100
plt.ylim(0,)  # Set y-axis to range from -0.35 to 0.35

plt.xlabel('Step')
plt.ylabel('Tracking Errors')
# plt.title('Simulation Results')
plt.show()
