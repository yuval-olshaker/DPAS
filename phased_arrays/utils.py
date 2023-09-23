import numpy as np
# Parameters
frequency = 1e9  # 1 GHz
c = 3e8  # Speed of light
wavelength = c / frequency # 0.3
k = 2 * np.pi / wavelength #wave number
epsilon = 1e-10
antenna_size = wavelength / 2


# Array dimensions
Nx = 10  # Number of elements along x
Ny = 10  # Number of elements along y
dx = 0.5 * wavelength  # Spacing along x
dy = 0.5 * wavelength  # Spacing along y
x = np.arange(Nx) * dx
y = np.arange(Ny) * dy

# Angles
THETA = np.linspace(-np.pi / 2, np.pi / 2, 400)
PHI = np.linspace(-np.pi, np.pi, 400)
azimuth_resolution = 1 # milliradian
azimuth_resolution = (azimuth_resolution / (1000 * 2 * np.pi)) * 360 # degrees
elevation_resolution = 2 # milliradian
elevation_resolution = (elevation_resolution / (1000 * 2 * np.pi)) * 360 # degrees
azimuth_range = 120 #degrees
elevation_range = 60 #degrees

# Pulse energy at each antenna
Pt_antenna = 1 # in watt per pulse
Pt_array = Pt_antenna * Nx * Ny

# Target
R = 100 * 1e3 #50 km
RCS = 10

# Noise
temperature = 290
kB = 1.38 * 1e-23
pulse_bandwidth = 50 * 1e6

# Params for SNR
G = 10000 # 40 DB
Ae = Nx * Ny * 1 # 1 meter spacing for antenna size. affective size is the whole size