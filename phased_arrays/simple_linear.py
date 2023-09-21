import numpy as np
import matplotlib.pyplot as plt

epsilon = 10**-9

def isotropic_pattern(theta):
    """Return the gain of an isotropic antenna (equal in all directions)."""
    return 1.0

# Define the gain function for a half-wave dipole antenna
def dipole_gain(theta):
    numerator = np.cos(np.pi/2 * np.cos(theta))
    denominator = np.sin(theta)
    return 1.643 * np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator!=0)**2 #the equation

# Gain function for slot antenna
def slot_gain(theta):
    return dipole_gain(np.pi/2 - theta)


def phased_array_pattern(theta, num_elements, d, wavelength, steering_angle):
    """Calculate the phased array pattern for a linear array."""
    n = np.arange(num_elements)
    delta_phi = (2 * np.pi / wavelength) * d * n * np.sin(steering_angle)
    array_factor = np.sum(np.exp(1j * (2 * np.pi / wavelength) * d * n[:, np.newaxis] * np.sin(theta) - 1j * delta_phi[:, np.newaxis]), axis=0)
    phased_array_gain = slot_gain(theta) * array_factor # the gain of the entire array
    phased_array_gain[phased_array_gain == 0] = epsilon # 0 value cannot have log
    return phased_array_gain

if __name__ == '__main__':
    # Parameters
    num_elements = 50
    d = 0.5  # Spacing (in terms of wavelength)
    wavelength = 1
    steering_angle = np.deg2rad(0)  # Desired steering angle in radians

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
    plt.ylim([-40, 40])  # Displaying from -40 dB to 0 dB
    plt.show()
