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
    print('total_power_per_scan_per_antenna: ' + str(total_power_per_scan_per_antenna(Pt_antenna, azimuth_resolution, elevation_resolution, azimuth_range, elevation_range) / 1000) + ' KW')
    print('SNR: ' + str(calculate_SNR(Pt_array, G, Ae, RCS, R, temperature, kB, pulse_bandwidth)) + ' DB')


