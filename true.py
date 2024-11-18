import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import odeint

# DC Motor Parameters
J = 0.01  # Moment of inertia
b = 0.1   # Motor friction
K = 0.01  # Motor constant
R = 1     # Electric resistance
L = 0.5   # Electric inductance

# State-space model for each agent
def state_space_model(x, t, u):
    A = np.array([[-b/J, K/J],
                  [-K/L, -R/L]])
    B = np.array([0, 1/L])
    dxdt = A @ x + B * u
    return dxdt

# MFAC control parameters
Kp = 1.0
Ki = 0.1
Kd = 0.01

# Desired positions
desired_position = np.array([1, 1])  # Desired positions for agent 1 and 2

# Initial states
x1 = np.array([0, 0])  # Initial state of agent 1
x2 = np.array([0, 0])  # Initial state of agent 2
u1 = 0
u2 = 0

# Time parameters
T = np.linspace(0, 10, 1000)

# Arrays to store results
x1_res = []
x2_res = []
u1_res = []
u2_res = []

# Control loop simulation
for t in T:
    # Current outputs
    y1 = x1[0]
    y2 = x2[0]

    # Errors
    e1 = desired_position[0] - y1
    e2 = desired_position[1] - y2

    # Control laws (MFAC)
    u1 = u1 + Kp * e1 + Ki * np.sum(e1) + Kd * (e1 - (desired_position[0] - y1))
    u2 = u2 + Kp * e2 + Ki * np.sum(e2) + Kd * (e2 - (desired_position[1] - y2))

    # Integrate the state-space model
    x1 = odeint(state_space_model, x1, [t, t+0.01], args=(u1,))[-1]
    x2 = odeint(state_space_model, x2, [t, t+0.01], args=(u2,))[-1]

    # Store results
    x1_res.append(x1)
    x2_res.append(x2)
    u1_res.append(u1)
    u2_res.append(u2)

# Convert results to numpy arrays
x1_res = np.array(x1_res)
x2_res = np.array(x2_res)
u1_res = np.array(u1_res)
u2_res = np.array(u2_res)

# Plot results
plt.figure(figsize=(14, 6))

plt.subplot(2, 2, 1)
plt.plot(T, x1_res[:, 0], label='Agent 1 Position')
plt.plot(T, np.ones_like(T) * desired_position[0], 'r--', label='Desired Position')
plt.xlabel('Time (s)')
plt.ylabel('Position')
plt.title('Agent 1 Position')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(T, u1_res, label='Agent 1 Control Input')
plt.xlabel('Time (s)')
plt.ylabel('Control Input')
plt.title('Agent 1 Control Input')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(T, x2_res[:, 0], label='Agent 2 Position')
plt.plot(T, np.ones_like(T) * desired_position[1], 'r--', label='Desired Position')
plt.xlabel('Time (s)')
plt.ylabel('Position')
plt.title('Agent 2 Position')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(T, u2_res, label='Agent 2 Control Input')
plt.xlabel('Time (s)')
plt.ylabel('Control Input')
plt.title('Agent 2 Control Input')
plt.legend()

plt.tight_layout()
plt.show()
