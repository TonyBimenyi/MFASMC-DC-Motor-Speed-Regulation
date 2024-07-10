import numpy as np

class MFA_controller:
    def __init__(self, rho=0.25 , eta=1, lambda_ = 1 , mu = 2.5):
        self.rho = rho
        self.eta = eta
        self.lambda_ = lambda_
        self.mu = mu

    def update_phi(self, k, phi, u, y):
        if k==1:
            return 2
        elif k == 2:
             return phi[k-1] + (self.eta * (u[k-1] - 0) / (self.mu + (u[k-1] - 0)**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - 0))
        else:
            return phi[k-1] + (self.eta * (u[k-1] - u[k-2]) / (self.mu + (u[k-1] - u[k-2])**2)) * (y[k] - y[k-1] - phi[k-1] * (u[k-1] - u[k-2]))
        
    def update_si(self, k, yd, y1, y2, y3, y4):
        if k == 1:
            return 0,0,0,0
        else: 
            return(yd[k] - 2*y1[k] + y4[k], 
                    y1[k] - 2*y2[k] + y3[k], 
                    y2[k] + yd[k] - 2*y3[k], 
                    y1[k] + y3[k] - 2*y4[k])
        
    def update_input(self, k, u, phi, si):
        if k == 1:
            return 0
        else:
            return u[k-1] + (self.rho * phi[k] / (self.lamda + phi[k]**2)) * si[k]