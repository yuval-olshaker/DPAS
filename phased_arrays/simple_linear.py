import numpy as np
import matplotlib.pyplot as plt

def isotropic_pattern(theta):
    """Return the gain of an isotropic antenna (equal in all directions)."""
    return 1.0

def phased_array_pattern(theta, num_elements, d, wavelength, steering_angle):
    """Calculate the phased array pattern for a linear array."""
    n = np.arange(num_elements)
    delta_phi = (2 * np.pi / wavelength) * d * n * np.sin(steering_angle)
    array_factor = np.sum(np.exp(1j * (2 * np.pi / wavelength) * d * n[:, np.newaxis] * np.sin(theta) - 1j * delta_phi[:, np.newaxis]), axis=0)
    return isotropic_pattern(theta) * array_factor

if __name__ == '__main__':
    # Parameters
    num_elements = 16
    d = 0.5  # Spacing (in terms of wavelength)
    wavelength = 1
    steering_angle = np.deg2rad(30)  # Desired steering angle in radians

    # Calculate patterns
    theta = np.linspace(-np.pi/2, np.pi/2, 1000)
    array_response = phased_array_pattern(theta, num_elements, d, wavelength, steering_angle)

    # Plot
    plt.figure(figsize=(10,6))
    plt.plot(np.rad2deg(theta), 20 * np.log10(np.abs(array_response)))
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Array Response (dB)')
    plt.title('Phased Array Antenna Simulation')
    plt.grid(True)
    plt.ylim([-40, 0])  # Displaying from -40 dB to 0 dB
    plt.show()
