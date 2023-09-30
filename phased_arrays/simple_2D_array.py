from Versatile_phased_array.utils import *


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
    return E

def dipole_simple_one_dim_pattern(theta):
    return np.cos(np.pi / 2 * np.cos(theta))

def isotopic_pattern(theta):
    return np.ones_like(theta)

def phased_array_pattern(x, y, kx, ky, phasex, phasey):
    """
    Compute the array factor for a 2D rectangular grid of antennas.
    """
    AF = np.zeros_like(kx, dtype=complex)

    for i in range(len(x)):
        for j in range(len(y)):
            AF += np.exp(1j * (kx * x[i] + ky * y[j] + phasex[i] + phasey[j]))

    return AF


def array_factor2(theta, phi, beta):
    af = np.zeros((theta.shape[0], phi.shape[0]))

    for i, t in enumerate(theta):
        for j, p in enumerate(phi):
            for m in range(Nx):
                for n in range(Ny):
                    af[i, j] += np.exp(
                        1j * (k * dx * m * np.sin(t) * np.cos(p) + k * dy * n * np.sin(t) * np.sin(p) - beta[m, n]))

    return af

def array_factor(theta, phi, beta):
    x_indices = np.arange(Nx)
    y_indices = np.arange(Ny)
    x_grid, y_grid, theta_grid, phi_grid = np.meshgrid(x_indices, y_indices, theta, phi, indexing='ij')

    phase = (k * dx * x_grid * np.sin(theta_grid) * np.cos(phi_grid) +
             k * dy * y_grid * np.sin(theta_grid) * np.sin(phi_grid) -
             beta[x_grid, y_grid])

    af = np.sum(np.exp(1j * phase), axis=(0, 1))

    return af

# The phases function
def calculate_beta(steer_theta, steer_phi):
    beta = np.zeros((Nx, Ny))
    for i in range(Nx):
        for j in range(Ny):
            beta[i, j] = -k * dx * i * np.sin(steer_theta) * np.cos(steer_phi) - k * dy * j * np.sin(steer_theta) * np.sin(steer_phi)

    return beta

def calculate_phases(theta, phi):
    phaseX = [2 * np.pi * dx * i * np.sin(phi) / wavelength for i in range(Nx)]
    phaseY = [2 * np.pi * dy * j * np.sin(theta) / wavelength for j in range(Ny)]
    return phaseX, phaseY

if __name__ == '__main__':
    # Phase distribution (for beam steering)
    steer_theta, steer_phi = np.radians(0), np.radians(0)  # Desired steering angles in radians
    beta = calculate_beta(steer_theta, steer_phi)

    # phasesX, phasesY = calculate_phases(steer_theta, steer_phi)
    # kx = k * np.outer(np.sin(THETA), np.cos(PHI))
    # ky = k * np.outer(np.sin(THETA), np.sin(PHI))
    # PA_pattern = phased_array_pattern(x, y, kx, ky, phasesX, phasesY)

    PA_pattern = array_factor(THETA, PHI, beta)
    AF = np.abs(PA_pattern)
    # AF2 = np.abs(array_factor2(THETA, PHI, beta))
    # print(np.array_equal(AF, AF2))
    # exit(0)
    IP = isotopic_pattern(THETA)
    DP = dipole_simple_one_dim_pattern(THETA)
    DP_gain = 1.643 * np.abs(DP) # the max gain is 1.643
    SA = IP
    total_pattern = AF * SA
    total_pattern[total_pattern == 0] = epsilon
    pattern_in_DB = 20 * np.log10(total_pattern)

    PHI_deg = np.rad2deg(PHI)
    THETA_deg = np.rad2deg(THETA)
    print('max gain DB: ' + str(np.max(pattern_in_DB)))
    a = np.unravel_index(np.argmax(pattern_in_DB), pattern_in_DB.shape)
    print('max db theta angle: ' + str(THETA_deg[a[0]]))
    print('max db phi angle: ' + str(PHI_deg[a[1]]))

    # Plot
    plt.figure()
    plt.pcolormesh(PHI_deg, THETA_deg, pattern_in_DB, shading='auto', vmin=-40, vmax=40)
    plt.colorbar(label='Radiation Pattern (dB)')
    plt.title('Planar Antenna Array with Dipole Antennas Radiation Pattern')
    plt.xlabel('Phi')
    plt.ylabel('Theta')
    plt.show()
