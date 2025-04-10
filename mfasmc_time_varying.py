import numpy as np
import matplotlib.pyplot as plt 

#define parameters 
d =5
rho = 1
eta = 1
lamda = 50
mu = 1
epsilon = 10**-5
alpha = 1
gamma1 = 0.15
gamma2 = 0.15
gamma3 = 0.45
gamma4 = 0.45

rt = 1024
m = 350
L = 200

yd = np.zeros(L + 1)
for k in range(L):
    yd[k] = 0.5 * np.sin(0.07 * np.pi * k) + 0.7 * np.cos(0.04 * np.pi * k)


plt.figure(figsize=(10, 6))  # Adjust figure size for better visibility

plt.plot(yd[:-1], '-y', linewidth=2.5, label=r'$y_d$')  

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