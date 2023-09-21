import matplotlib.pyplot as plt
from phased_arrays.utils import *


def dipole_pattern(theta, phi, L, k):
    """"
    Compute the radiation pattern for a dipole antenna oriented along the x-axis.
    """
    # Radiation pattern for a dipole of length L oriented along the x-axis
    E = np.cos(phi) * (np.cos(k * L * np.cos(theta) / 2) - np.cos(k * L / 2)) / np.sin(theta)

    # Handle the singularity at theta = 0 and theta = pi
    E[np.isnan(E)] = 0
    return E


def phased_array_pattern(x, y, kx, ky, phasex, phasey):
    """
    Compute the array factor for a 2D rectangular grid of antennas.
    """
    Nx, Ny = phasex.shape
    AF = np.zeros_like(kx, dtype=complex)

    for i in range(Nx):
        for j in range(Ny):
            AF += np.exp(1j * (kx * x[i] + ky * y[j] + phasex[i , j] + phasey[i, j]))

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
    phasex = np.array([
        [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
        [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
        [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
        [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]
    ])

    phasey = np.array([
        [0, -np.pi / 4, -np.pi / 2, -3 * np.pi / 4],
        [0, -np.pi / 4, -np.pi / 2, -3 * np.pi / 4],
        [0, -np.pi / 4, -np.pi / 2, -3 * np.pi / 4],
        [0, -np.pi / 4, -np.pi / 2, -3 * np.pi / 4]
    ])
    # # Linear phase progression in x-direction (e.g., 45 degrees per element)
    # phasex = np.deg2rad(45) * np.arange(Nx)
    # # Quadratic phase progression in y-direction (e.g., 30 degrees per element squared)
    # phasey = np.deg2rad(30) * (np.arange(Ny) ** 2)
    print(phasex)
    # Compute radiation pattern
    theta = np.linspace(0, np.pi, 400)
    phi = np.linspace(0, 2 * np.pi, 400)
    THETA, PHI = np.meshgrid(theta, phi)

    kx = k * np.sin(THETA) * np.cos(PHI)
    ky = k * np.sin(THETA) * np.sin(PHI)

    pattern = phased_array_pattern(x, y, kx, ky, phasex, phasey)
    AF = np.abs(pattern)
    DP = dipole_pattern(THETA, PHI, antenna_size, k)

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
