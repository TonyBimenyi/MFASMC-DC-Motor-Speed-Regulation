import numpy as np
from Controllers.MFA_controller import MFA_controller
from Controllers.SM_controller import SM_controller





class MFASM_controller:
    def __init__(self, gamma=20, eta=1, lambda_ = 1 , mu = 2.5):
        self.gamma = gamma
        self.MFA =  MFA_controller()
        self.SM = SM_controller()
        self.eta =eta
        self.lambda_ = lambda_
        self.mu = mu

    def update_phi(self, k, phi, u, y):
        if k==1:
            return 2
        elif k == 2:
             return phi[k-1] + (self.eta * (u[k-1] - 0) / (self.mu + (u[k-1] - 0)**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - 0))
        else:
            return phi[k-1] + (self.eta * (u[k-1] - u[k-2]) / (self.mu + (u[k-1] - u[k-2])**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - u[k-2]))

    def update_u(self,k):
        if k == 1:
            return 0
        else: 
            return self.MFA +  self.SM