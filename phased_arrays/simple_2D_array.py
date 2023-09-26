from phased_arrays.utils import *


def paper_gain(theta,b, k):
    term = k * theta * b / 2
    gain = (k * b) ** 2 / np.pi * (np.sin(term) ** 2) / (term ** 2)

    return gain

def dipole_pattern(theta, phi, L, k):
    """"
    Compute the radiation pattern for a dipole antenna oriented along the x-axis.
    """
    # Radiation pattern for a dipole of length L oriented along the x-axis
    E = np.cos(phi) * ((np.cos(k * L * np.cos(theta) / 2) - np.cos(k * L / 2)) / np.sin(theta))

    # Handle the singularity at theta = 0 and theta = pi
    E[np.isnan(E)] = 0
    return np.abs(E)


def phased_array_pattern(x, y, kx, ky, phasex, phasey):
    """
    Compute the array factor for a 2D rectangular grid of antennas.
    """
    AF = np.zeros_like(kx, dtype=complex)

    for i in range(len(x)):
        for j in range(len(y)):
            AF += np.exp(1j * (kx * x[i] + ky * y[j] + phasex[i] + phasey[j]))

    return AF

def calculate_phases(Nx, Ny, dx, dy, lambda_, theta_x, theta_y):
    phaseX = [2 * np.pi * dx * i * np.sin(theta_x) / lambda_ for i in range(Nx)]
    phaseY = [2 * np.pi * dy * j * np.sin(theta_y) / lambda_ for j in range(Ny)]
    return phaseX, phaseY

if __name__ == '__main__':
    # Phase distribution (for beam steering)
    theta, phi = np.radians(0), np.radians(0)  # Desired steering angles in radians
    phasesX, phasesY = calculate_phases(Nx, Ny, dx, dy, wavelength, theta, phi)

    kx = k * np.outer(np.sin(THETA), np.cos(PHI))
    ky = k * np.outer(np.sin(THETA), np.sin(PHI))

    PA_pattern = phased_array_pattern(x, y, kx, ky, phasesX, phasesY)
    AF = np.abs(PA_pattern)
    DP = dipole_pattern(THETA, PHI, antenna_size, k)
    DP_gain = 1.643 * DP # the max gain is 1.643

    total_pattern = AF * DP_gain
    total_pattern[total_pattern == 0] = epsilon
    pattern_in_DB = 20 * np.log10(total_pattern)
    print(np.max(pattern_in_DB))
    print(np.argmax(pattern_in_DB))

    # Plot
    plt.figure()
    PHI_deg = np.rad2deg(PHI)
    THETA_deg = np.rad2deg(THETA)
    plt.pcolormesh(PHI_deg, THETA_deg, pattern_in_DB, shading='auto', vmin=-40, vmax=40)
    plt.colorbar(label='Radiation Pattern (dB)')
    plt.title('Planar Antenna Array with Dipole Antennas Radiation Pattern')
    plt.xlabel('Phi')
    plt.ylabel('Theta')
    plt.show()
