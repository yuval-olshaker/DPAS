from phased_arrays.utils import *

def radar_equation(Pt, G, Ae, sigma, R):
    """
    Calculate the received power using the radar equation.

    Parameters:
    - Pt: Transmitted power (W)
    - G: Antenna gain (dimensionless, not in dB)
    - Ae: Effective aperture (m^2)
    - sigma: Radar cross-section (m^2)
    - R: Range to the target (m)

    Returns:
    - Pr: Received power (W)
    """
    Pr = (Pt * G * Ae * sigma) / ((4 * np.pi)**2 * R**4)
    return Pr

def noise_power(temperature, kB, pulse_bandwidth):
    return temperature * kB * pulse_bandwidth

def calculate_SNR(Pt_array, G, Ae, RCS, R, temperature, kB, pulse_bandwidth):
    pure_SNR = radar_equation(Pt_array, G, Ae, RCS, R) / noise_power(temperature, kB, pulse_bandwidth)
    return 10 * np.log10(pure_SNR)

def total_power_per_scan_per_antenna(Pt_antenna, azimuth_resolution, elevation_resolution, azimuth_range, elevation_range):
    return Pt_antenna * (azimuth_range / azimuth_resolution) * (elevation_range / elevation_resolution)

if __name__ == '__main__':
    # change antenna Pt and drones number
    points_num = 20
    Pt_antenna = np.linspace(0.001, 0.1, points_num)
    Nx = np.linspace(5, 100, points_num)
    Ny = Nx
    SNRs = np.zeros((points_num, points_num))
    for i in range(points_num):
        for j in range(points_num):
            Pt_array = Pt_antenna[i] * Nx[j] * Ny[j]
            Ae = Nx[j] * Ny[j] * 1
            SNRs[i, j] = calculate_SNR(Pt_array, G, Ae, RCS, R, temperature, kB, pulse_bandwidth)
    print(SNRs)

    # print('total_power_per_scan_per_antenna: ' + str(total_power_per_scan_per_antenna(Pt_antenna, azimuth_resolution, elevation_resolution, azimuth_range, elevation_range) / 1000) + ' KW')
    # print('SNR: ' + str(calculate_SNR(Pt_array, G, Ae, RCS, R, temperature, kB, pulse_bandwidth)) + ' DB')


