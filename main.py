import numpy as np

from Controllers.MFASM import MFASM_controller
from Controllers.MFA_controller import MFA_controller
from Reference.trajectory import generate_trajectory
from Output.output import Output
from Estimator.estimator import initialize_phi
from Plot.plotter import plot_results

controller = MFASM_controller()
mfa_ = MFA_controller()
phi = initialize_phi()
yd = generate_trajectory(type='cosine_squared')
u = [np.zeros((400, 1)) for _ in range(4)]
y = [np.zeros((401, 1)) for _ in range(4)]
si = [np.zeros((400, 1)) for _ in range(4)]

for k in range(400):
    # Update phi values
    for i in range(4):
        phi[i][k] = mfa_.update_phi(k, phi[i], u[i], y[i])

    # Update si values
    si1, si2, si3, si4 = mfa_.update_si(k, yd, y[0], y[1], y[2], y[3])
    si[0][k], si[1][k], si[2][k], si[3][k] = si1, si2, si3, si4
    
    # Update input values
    for i in range(4):
        u[i][k] = controller.update_u(k,phi,si,u,yd,y)

    # Update y values
    y1_next, y2_next, y3_next, y4_next = Output(k, y)
    y[0][k+1], y[1][k+1], y[2][k+1], y[3][k+1] = y1_next, y2_next, y3_next, y4_next


plot_results(yd)