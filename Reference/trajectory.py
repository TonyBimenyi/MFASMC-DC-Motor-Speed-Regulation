import numpy as np

def generate_trajectory():
    yd = np.zeros(400)
    yd = 50 * np.cos(200 * np.pi /400) + 50

    return yd