import numpy as np

class SM_controller:
    def __ini__(self, epsilon = 10**-4 , alpha = 0.5, T=3):
        self.epsilon = epsilon
        self.alpha = alpha
        self.T = T

    def update_phi(self, k, phi, u, y):
        if k==1:
            return 2
        elif k == 2:
             return phi[k-1] + (self.eta * (u[k-1] - 0) / (self.mu + (u[k-1] - 0)**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - 0))
        else:
            return phi[k-1] + (self.eta * (u[k-1] - u[k-2]) / (self.mu + (u[k-1] - u[k-2])**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - u[k-2]))
        
    def update_input(self, k, u, phi, yd, y):
        if k ==1:
            return 0
        else:
            return u[k-1] + (yd[k+1] - y[k] + self.alpha * yd[k+1] - y[k] + self.epsilon * self.T * np.sign(k))