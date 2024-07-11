import numpy as np
from scipy import signal

def generate_trajectory(type='step', t=np.linspace(0, 1, 100)):
    if type == 'step':
        yd = np.zeros(len(t))
        yd[:len(t)//2] = 1
        yd[len(t)//2:] = 2
    elif type == 'sinusoidal':
        yd = np.sin(t)
    elif type == 'ramp':
        yd = np.linspace(0, 10, len(t))
    elif type == 'quadratic':
        yd = t**2
    elif type == 'exponential':
        yd = np.exp(t)
    elif type == 'sine_with_amplitude_modulation':
        yd = np.sin(t) * np.sin(0.2 * t)
    elif type == 'cubic':
        yd = t**3
    elif type == 'sawtooth':
        yd = signal.sawtooth(t)
    elif type == 'triangular':
        yd = signal.sawtooth(t, 0.5)
    elif type == 'custom_piecewise_linear':
        yd = np.zeros(len(t))
        yd[:len(t)//4] = np.linspace(0, 1, len(t)//4)
        yd[len(t)//4:len(t)//2] = np.linspace(1, 3, len(t)//4)
        yd[len(t)//2:3*len(t)//4] = np.linspace(3, 0, len(t)//4)
        yd[3*len(t)//4:] = np.linspace(0, 2, len(t)//4)
    elif type == 'square_wave':
        yd = signal.square(t)
    elif type == 'sine_squared':
        yd = np.sin(t)**2
    elif type == 'cosine':
        yd = np.cos(t)
    elif type == 'cosine_squared':
        yd = np.cos(t)**2
    elif type == 'logarithmic':
        yd = np.log(t + 1)  # log(0) is undefined, hence t + 1
    elif type == 'inverse':
        yd = 1 / (t + 1)  # 1/0 is undefined, hence t + 1
    elif type == 'heaviside_step':
        yd = np.heaviside(t - 2, 1)  # Step function that steps at t = 2
    elif type == 'sinusoidal_with_decay':
        yd = np.sin(t) * np.exp(-0.5 * t)
    elif type == 'cosine_with_decay':
        yd = np.cos(t) * np.exp(-0.5 * t)
    elif type == 'bouncing_ball':
        yd = np.abs(np.sin(t))
    elif type == 'sigmoid':
        yd = 1 / (1 + np.exp(-t))
    elif type == 'linear_piecewise':
        yd = np.piecewise(t, [t < 1, (t >= 1) & (t < 2), (t >= 2) & (t < 3), t >= 3], [lambda t: t, lambda t: 2 - t, lambda t: t - 2, lambda t: 4 - t])
    elif type == 'parabolic':
        yd = t**2 - 4*t + 4
    elif type == 'sinusoidal_half_wave_rectified':
        yd = np.maximum(0, np.sin(t))
    elif type == 'cosine_half_wave_rectified':
        yd = np.maximum(0, np.cos(t))
    elif type == 'full_wave_rectified_sine':
        yd = np.abs(np.sin(t))
    elif type == 'full_wave_rectified_cosine':
        yd = np.abs(np.cos(t))
    elif type == 'linear_ramp_down':
        yd = -np.linspace(0, 10, len(t))
    elif type == 'sinusoidal_with_random_noise':
        yd = np.sin(t) + np.random.normal(0, 0.1, len(t))
    elif type == 'cosine_with_random_noise':
        yd = np.cos(t) + np.random.normal(0, 0.1, len(t))
    else:
        raise ValueError("Unsupported trajectory type")
    return yd
    