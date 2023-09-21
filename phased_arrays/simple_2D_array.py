import matplotlib.pyplot as plt
from phased_arrays.utils import *


def dipole_pattern(theta):
    """
    Compute the radiation pattern for a half-wave dipole antenna.
    """
    # Avoid division by zero at theta = 0 and theta = pi
    theta = np.clip(theta, epsilon, np.pi - epsilon)
    return np.cos(np.pi / 2 * np.cos(theta)) / np.sin(theta)


def phased_array_pattern(x, y, kx, ky, phasex, phasey):
    """
    Compute the array factor for a 2D rectangular grid of antennas.
    """
    Nx = len(x)
    Ny = len(y)

    AF = 0
    for i in range(Nx):
        for j in range(Ny):
            AF += np.exp(1j * (kx * x[i] + ky * y[j] + phasex[i] + phasey[j]))

    return AF

if __name__ == '__main__':
    # Array dimensions
    Nx = 5  # Number of elements along x
    Ny = 5  # Number of elements along y
    dx = 0.5 * wavelength  # Spacing along x
    dy = 0.5 * wavelength  # Spacing along y

    x = np.arange(Nx) * dx
    y = np.arange(Ny) * dy

    # Phase distribution (for beam steering)
    # Linear phase progression in x-direction (e.g., 45 degrees per element)
    phasex = np.deg2rad(45) * np.arange(Nx)
    # Quadratic phase progression in y-direction (e.g., 30 degrees per element squared)
    phasey = np.deg2rad(30) * (np.arange(Ny) ** 2)
    print(phasex)
    # Compute radiation pattern
    theta = np.linspace(0, np.pi, 400)
    phi = np.linspace(0, 2 * np.pi, 400)
    THETA, PHI = np.meshgrid(theta, phi)

    kx = k * np.sin(THETA) * np.cos(PHI)
    ky = k * np.sin(THETA) * np.sin(PHI)

    pattern = phased_array_pattern(x, y, kx, ky, phasex, phasey)
    AF = np.abs(pattern)
    DP = dipole_pattern(THETA)

    total_pattern = AF * DP

    # Plot
    plt.figure()
    pattern_in_DB = 20 * np.log10(total_pattern)
    plt.pcolormesh(PHI, THETA, pattern_in_DB, shading='auto', vmin=-40, vmax=40)
    plt.colorbar(label='Radiation Pattern (dB)')
    plt.title('Planar Antenna Array with Dipole Antennas Radiation Pattern')
    plt.xlabel('Phi (radians)')
    plt.ylabel('Theta (radians)')
    plt.show()
