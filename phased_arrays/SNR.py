from phased_arrays.utils import *
import pandas as pd

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

def noise_power():
    return temperature * kB * pulse_bandwidth

def calculate_SNR(Pt_array, G, Ae, RCS, R):
    pure_SNR = radar_equation(Pt_array, G, Ae, RCS, R) / noise_power()
    return 10 * np.log10(pure_SNR)

def total_power_per_scan_per_antenna(Pt_antenna, azimuth_resolution, elevation_resolution, azimuth_range, elevation_range):
    return Pt_antenna * (azimuth_range / azimuth_resolution) * (elevation_range / elevation_resolution)

def calculate_SNRs_for_changing_drones_num_and_Pt():
    # change antenna Pt and drones number
    Pt_antenna = np.linspace(1, 20, points_num)
    Nx = np.linspace(4, 20, points_num)
    Nx = np.round(Nx)
    Ny = Nx
    SNRs = np.zeros((points_num, points_num))
    for i in range(points_num):
        for j in range(points_num):
            Pt_array = Pt_antenna[i] * Nx[j] * Ny[j]
            Ae = Nx[j] * Ny[j] * 1
            SNRs[i, j] = calculate_SNR(Pt_array, G, Ae, RCS, R)

    return SNRs, Nx, Pt_antenna

def plot_table(SNRs, Nx, Pt_antenna):
    # Create a 2D array
    data = np.around(SNRs, 2)

    # Create column and row labels
    columns = np.around(Nx, 2).tolist()
    rows = np.around(Pt_antenna, 3)

    # Parameters to display on the side
    params = {
        'Parameter 1': 'Value 1',
        'Parameter 2': 'Value 2',
        'Parameter 3': 'Value 3'
    }

    fig, ax = plt.subplots(figsize=(8, 6))

    # Hide axes
    ax.axis('off')

    # Add the table
    table = ax.table(cellText=data, colLabels=columns, rowLabels=rows, cellLoc='center', loc='center')

    # Adjust table position to fit within the figure
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Add title
    plt.title("SNR as function of drones number and transmitted power")

    # Display parameters on the side
    for i, (key, value) in enumerate(params.items()):
        ax.text(1.05, 0.9 - i * 0.1, f"{key}: {value}", transform=ax.transAxes)
    # Add titles for columns and rows
    ax.text(0.5, 0.04, "Drones Num per Axis", size=12, ha="center", transform=fig.transFigure)
    ax.text(0.02, 0.5, "Pt per Drone Pulse (W)", size=12, va="center", rotation="vertical", transform=fig.transFigure)

    plt.tight_layout()
    plt.show()

def to_exel(SNRs, Nx, Pt_antenna):
    matrix = np.zeros((points_num + 2, points_num + 1))
    for i in range(points_num):
        matrix[0, i + 1] = Nx[i]
        matrix[1, i + 1] = Nx[i] ** 2
        matrix[i + 2, 0] = Pt_antenna[i]
        for j in range(points_num):
            matrix[i + 2, j + 1] = SNRs[i, j]

    matrix = np.around(matrix, 2) # 2 digits after the point
    df = pd.DataFrame(matrix)
    df.to_excel('output.xlsx', index=False, header=False)


if __name__ == '__main__':
    points_num = 10
    SNRs, Nx, Pt_antenna = calculate_SNRs_for_changing_drones_num_and_Pt()
    # Export to excel table
    to_exel(SNRs, Nx, Pt_antenna)
    # Plot table
    # plot_table(SNRs, Nx, Pt_antenna)

    # # Plot Fig
    # plt.figure()
    # plt.pcolormesh(Nx, Pt_antenna, SNRs, shading='auto', vmin=-25, vmax=25)
    # plt.colorbar(label='SNR (dB)')
    # plt.title('SNR as function of drones number and transmitted power')
    # plt.xlabel('Drones Num per Axis')
    # plt.ylabel('Pt per Drone Pulse (W)')
    # plt.show()
    # print(np.max(SNRs))
    # print(total_power_per_scan)
    # print((azimuth_range / azimuth_resolution) * (elevation_range / elevation_resolution))
    # print('total_power_per_scan_per_antenna: ' + str(total_power_per_scan_per_antenna(Pt_antenna, azimuth_resolution, elevation_resolution, azimuth_range, elevation_range) / 1000) + ' KW')
    # print('SNR: ' + str(calculate_SNR(Pt_array, G, Ae, RCS, R, temperature, kB, pulse_bandwidth)) + ' DB')


