import numpy as np
# Parameters
frequency = 1e9  # 3 GHz
c = 3e8  # Speed of light
wavelength = c / frequency
k = 2 * np.pi / wavelength
epsilon = 1e-10