import numpy as np
import matplotlib.pyplot as plt 

d = 5
rho = 0.5
eta = 1.2
lamda = 7
mu = 10
epsilon = 10**-5
alpha = 70
T = 0.1
gamma1 = 0.15
gamma2 = 0.45
gamma3 = 0.15
gamma4 = 0.45
alpha = 100
omega = 0.2
sigma = 600
tau = 0.3

rT = 1024
m = 1024
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

xi1 = np.zeros((L, 1))
xi2 = np.zeros((L, 1))
xi3 = np.zeros((L, 1))
xi4 = np.zeros((L, 1))

s1 = np.zeros((L, 1))
s2 = np.zeros((L, 1))
s3 = np.zeros((L, 1))
s4 = np.zeros((L, 1))

# Simulation loop
for k in range(1, L-1):
    if k == 1:
        phi1[1] = 1.5
        phi2[1] = 1.5
        phi3[1] = 1.5
        phi4[1] = 1.5
    elif k == 2:
        phi1[k] = phi1[k - 1] + eta * u1[k - 1] / (mu + u1[k - 1]**2) * (y1[k] - y1[k - 1] - phi1[k - 1] * u1[k - 1])
        phi2[k] = phi2[k - 1] + eta * u2[k - 1] / (mu + u2[k - 1]**2) * (y2[k] - y2[k - 1] - phi2[k - 1] * u2[k - 1])
        phi3[k] = phi3[k - 1] + eta * u3[k - 1] / (mu + u3[k - 1]**2) * (y3[k] - y3[k - 1] - phi3[k - 1] * u3[k - 1])
        phi4[k] = phi4[k - 1] + eta * u4[k - 1] / (mu + u4[k - 1]**2) * (y4[k] - y4[k - 1] - phi4[k - 1] * u4[k - 1])
    else:
        phi1[k] = phi1[k - 1] + (eta * (u1[k - 1] - u1[k - 2]) / (mu + ((u1[k - 1] - u1[k - 2]))**2)) * (y1[k] - y1[k - 1] - phi1[k - 1] * (u1[k - 1] - u1[k - 2]))
        phi2[k] = phi2[k - 1] + (eta * (u2[k - 1] - u2[k - 2]) / (mu + ((u2[k - 1] - u2[k - 2]))**2)) * (y2[k] - y2[k - 1] - phi2[k - 1] * (u2[k - 1] - u2[k - 2]))
        phi3[k] = phi3[k - 1] + (eta * (u3[k - 1] - u3[k - 2]) / (mu + ((u3[k - 1] - u3[k - 2]))**2)) * (y3[k] - y3[k - 1] - phi3[k - 1] * (u3[k - 1] - u3[k - 2]))
        phi4[k] = phi4[k - 1] + (eta * (u4[k - 1] - u4[k - 2]) / (mu + ((u4[k - 1] - u4[k - 2]))**2)) * (y4[k] - y4[k - 1] - phi4[k - 1] * (u4[k - 1] - u4[k - 2]))

    # Stability checks moved here
    if k > 1 and (abs(phi1[k]) <= epsilon or abs(u1[k - 1] - u1[k - 2]) <= epsilon or np.sign(phi1[k]) != np.sign(phi1[1])):
        phi1[k] = phi1[1]
    
    if k > 1 and (abs(phi2[k]) <= epsilon or abs(u2[k - 1] - u2[k - 2]) <= epsilon or np.sign(phi2[k]) != np.sign(phi2[1])):
        phi2[k] = phi2[1]

    if k > 1 and (abs(phi3[k]) <= epsilon or abs(u3[k - 1] - u3[k - 2]) <= epsilon or np.sign(phi3[k]) != np.sign(phi3[1])):
        phi3[k] = phi3[1]

    if k > 1 and (abs(phi4[k]) <= epsilon or abs(u4[k - 1] - u4[k - 2]) <= epsilon or np.sign(phi4[k]) != np.sign(phi4[1])):
        phi4[k] = phi4[1]

   

    xi1[k] = yd[k] - 2 * y1[k] + y4[k]
    xi2[k] = y1[k] - 2 * y2[k] + y3[k]
    xi3[k] = y2[k] + yd[k] - 2 * y3[k]
    xi4[k] = y1[k] + y3[k] - 2 * y4[k]
    # xi1[k] = yd[k]-y1[k]
    # xi2[k] = yd[k]-y2[k]
    # xi3[k] = yd[k]-y3[k]
    # xi4[k] = yd[k]-y4[k]

    # xi1[k+1] = yd[k+1] - 2 * y1[k+1] + y4[k+1]
    # xi2[k+1] = y1[k+1] - 2 * y2[k+1] + y3[k+1]
    # xi3[k+1] = y2[k+1] + yd[k+1] - 2 * y3[k+1]
    # xi4[k+1] = y1[k+1] + y3[k+1] - 2 * y4[k+1]

    
    if k == 0:
        s1[0] = 0
        s2[0] = 0
        s3[0] = 0
        s4[0] = 0
    else:
        s1[k] = alpha * xi1[k] - xi1[k-1]
        s2[k] = alpha * xi2[k] - xi2[k-1]
        s3[k] = alpha * xi3[k] - xi3[k-1]
        s4[k] = alpha * xi4[k] - xi4[k-1]


    
    


    if k == 1:
        mfa1[0] = 0
        mfa2[0] = 0
        mfa3[0] = 0
        mfa4[0] = 0
    else:
        mfa1[k] = mfa1[k - 1] + (rho * phi1[k]) / (lamda + (phi1[k])**2) * xi1[k]
        mfa2[k] = mfa2[k - 1] + (rho * phi2[k]) / (lamda + (phi2[k])**2) * xi2[k]
        mfa3[k] = mfa3[k - 1] + (rho * phi3[k]) / (lamda + (phi3[k])**2) * xi3[k]
        mfa4[k] = mfa4[k - 1] + (rho * phi4[k]) / (lamda + (phi4[k])**2) * xi4[k]

    if k == 1:
        sm1[0] = 0
        sm2[0] = 0
        sm3[0] = 0
        sm4[0] = 0
    else:
        # sm1[k] = sm1[k-1] + ((omega * phi1[k])/(xigma+abs(phi1[k])**2) * ((alpha*(xi1[k])-xi1[k])/alpha*(y4[k+1]+yd[k+1]))-y1[k]+tau*np.sign(s1[k]))

        # sm2[k] = sm2[k-1] + ((omega * phi2[k])/(sigma+abs(phi2[k])**2) * ((alpha*(xi2[k])-xi2[k])/alpha*(y1[k+1]+y3[k+1]+0))-y2[k]+tau*np.sign(s2[k]))

        # sm3[k] = sm3[k-1] + ((omega * phi3[k])/(sigma+abs(phi3[k])**2) * ((alpha*(xi3[k])-xi3[k])/alpha*(y2[k+1]+yd[k+1]))-y3[k]+tau*np.sign(s3[k]))

        # sm4[k] = sm4[k-1] + ((omega * phi4[k])/(sigma+abs(phi4[k])**2) * ((alpha*(xi4[k])-xi4[k])/alpha*(y1[k+1]+y3[k+1]+0))-y4[k]+tau*np.sign(s4[k]))

        

        sm1[k] = sm1[k-1] + (omega * phi1[k]) / (sigma + (phi1[k])**2) * ((xi1[k]+(y4[k]-y4[k-1])+(yd[k+1]-yd[k]))/(1+1)-(xi1[k])/(alpha*(1+1))+tau * np.sign(s1[k]))

        sm2[k] = sm2[k-1] + (omega * phi2[k]) / (sigma + (phi2[k])**2) * ((xi2[k]+(y1[k]-y1[k-1])+(y3[k]-y3[k-1])+(yd[k+1]-yd[k]))/(2)-(xi2[k])/(alpha*(2))+tau * np.sign(s2[k]))


        sm3[k] = sm3[k-1] + (omega * phi3[k]) / (sigma + (phi3[k])**2) * ((xi3[k]+(y2[k]-y2[k-1])+(yd[k+1]-yd[k]))/(1+1)-(xi3[k])/(alpha*(1+1))+tau * np.sign(s3[k]))

        sm4[k] = sm4[k-1] + (omega * phi4[k]) / (sigma + (phi4[k])**2) * ((xi4[k]+(y1[k]-y1[k-1])+(y3[k]-y3[k-1])+(yd[k+1]-yd[k]))/(1+1)-(xi4[k])/(alpha*(1+1))+tau * np.sign(s4[k]))
       

    
    
        
    if k == 1:
        u1[0] = 0
        u2[0] = 0
        u3[0] = 0
        u4[0] = 0
    else:
        u1[k] = mfa1[k] + gamma1 * sm1[k]
        u2[k] = mfa2[k] + gamma2 * sm2[k]
        u3[k] = mfa3[k] + gamma3 * sm3[k]
        u4[k] = mfa4[k] + gamma4 * sm4[k]

    y1[0] = 0.1
    y2[0] = 0.1
    y3[0] = 0.1
    y4[0] = 0.1
    # y1(k + 1) = 3*y1(k) * u1(k) / (15 + y1(k)**2) + 1.1*u1(k)
    # y2(k + 1) = 5*y2(k) * u2(k) / (355 + y2(k)**2) + 1.3*u2(k)
    # y3(k + 1) = 3*y3(k) * u3(k) / (15 + y3(k)**2) + 1.1*u3(k)
    # y4(k + 1) = 5*y4(k) * u4(k) / (355 + y4(k)**2) + 1.3*u4(k)

    # y1[k+1] = (y1[k] * u1[k])/(1+y1[k]**2) + 0.5*u1[k]
    # y2[k+1] = (y2[k] * u2[k])/(1+y2[k]**2) + 0.9*u2[k]
    # y3[k+1] = (y3[k] * u3[k])/(1+y3[k]**2) + 0.5*u3[k]
    # y4[k+1] = (y4[k] * u4[k])/(1+y4[k]**2) + 0.9*u4[k]

    y1[k+1] = (m / (rT * 0.3)) * u1[k]
    y2[k+1] = (m / (rT * 0.3)) * u2[k]
    y3[k+1] = (m / (rT * 0.3)) * u3[k]
    y4[k+1] = (m / (rT * 0.3)) * u4[k]


    

    
        

plt.figure(figsize=(10, 6))  # Adjust figure size for better visibility

# plt.plot(yd[:-1], '-y', linewidth=2, label=r'$y_d$')  
# plt.plot(y1[:-1], '--r', linewidth=2, label=r'$y_1$')
# plt.plot(y2[:-1], '-.b', linewidth=2, label=r'$y_2$')
# plt.plot(y3[:-1], '--k', linewidth=2, label=r'$y_3$')
# plt.plot(y4[:-1], '-.g', linewidth=2, label=r'$y_4$')
plt.plot(yd[:-1], '-y', label=r'$y_d$')  
plt.plot(y1[:-1], '--r', label=r'$y_1$')
plt.plot(y2[:-1], '-.b', label=r'$y_2$')
plt.plot(y3[:-1], '--k', label=r'$y_3$')
plt.plot(y4[:-1], '-.g', label=r'$y_4$')

# Enlarge font size for axis labels and legend
plt.xlabel('Time step', fontsize=14)
plt.ylabel('Output', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,L)  # Set y-axis to range from -0.35 to 0.35

# Improve legend readability
plt.legend(fontsize=14, loc='best')
# plt.title('Tracking performance',fontsize=15, fontweight='bold',fontname=font_style)

plt.grid(False)  # Add grid for better readability
plt.tight_layout()  # Adjust layout for clarity

plt.show()