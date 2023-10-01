import numpy as np
import matplotlib.pyplot as plt

# Parameters
frequency = 1e9  # 1 GHz
c = 3e8  # Speed of light
wavelength = c / frequency # 0.3
k = 2 * np.pi / wavelength #wave number
epsilon = 1e-10
antenna_size = wavelength / 2


# Array dimensions
Nx = 10  # Number of elements along x
Ny = Nx  # Number of elements along y
dx = 0.5 * wavelength  # Spacing along x
dy = 0.5 * wavelength  # Spacing along y
x = np.arange(Nx) * dx
y = np.arange(Ny) * dy

# Angles
azimuth_range = np.radians(240) # 120 degrees
elevation_range = np.radians(120) # 60 degrees
PHI = np.linspace(- azimuth_range / 2, azimuth_range / 2, 50)
THETA = np.linspace(-elevation_range / 2, elevation_range / 2, 50)
azimuth_resolution = 1 # milliradian
azimuth_resolution = (azimuth_resolution / (1000 * 2 * np.pi)) * 360 # degrees
elevation_resolution = 2 # milliradian
elevation_resolution = (elevation_resolution / (1000 * 2 * np.pi)) * 360 # degrees

# Pulse energy at each antenna
Pt_antenna = 1 # in watt per pulse
Pt_array = Pt_antenna * Nx * Ny

# Target
R = 100 * 1e3
RCS = 10

# Noise
temperature = 290
kB = 1.38 * 1e-23
pulse_bandwidth = 50 * 1e6

# Params for SNR
G = 5000 # 37 DB
Ae = Nx * Ny * 1 # 1 meter spacing for antenna size. affective size is the whole size

# Errors
pos_error_sigma_per_axis = 1 / 100 # 0.5 cm error per axis

def distance(pos1, pos2):
    return np.linalg.norm(np.array(pos1) - np.array(pos2))

def normalize_angle(angle):
    """
    Normalize an angle to the range [-pi, pi].

    Parameters:
    - angle (float): The angle in radians.

    Returns:
    - float: The normalized angle in the range [-pi, pi].
    """
    normalized_angle = angle % (2 * np.pi)
    if normalized_angle > np.pi:
        normalized_angle -= 2 * np.pi
    return normalized_angle