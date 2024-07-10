import numpy as np

def initialize_phi():
    phi = [np.zeros((400, 1)) for _ in range(4)]
    return phi