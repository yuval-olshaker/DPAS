import numpy as np
# Parameters
frequency = 1e9  # 1 GHz
c = 3e8  # Speed of light
wavelength = c / frequency # 0.3
k = 2 * np.pi / wavelength #wave number
epsilon = 1e-10
antenna_size = wavelength / 2


# Array dimensions
Nx = 5  # Number of elements along x
Ny = 5  # Number of elements along y
dx = 0.5 * wavelength  # Spacing along x
dy = 0.5 * wavelength  # Spacing along y
x = np.arange(Nx) * dx
y = np.arange(Ny) * dy

# Angles
theta = np.linspace(0, np.pi, 400)
phi = np.linspace(0, 2 * np.pi, 400)
THETA, PHI = np.meshgrid(theta, phi)

# Pulse energy at each antenna
Pt_antenna = 40 #in watt
Pt_array = Pt_antenna * Nx * Ny

# Target
R = 50 * 10e3 #50 km
RCS = 10