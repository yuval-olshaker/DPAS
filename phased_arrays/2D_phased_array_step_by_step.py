from phased_arrays.utils import *


def calculate_radiation_at_target(antennas, target_pos):
    """
     Calculate the radiation power densities and phases at a target position from a list of antennas.

     Parameters:
     - antennas (list): A list of antenna objects. Each antenna object should have a method
                        `radiation_at_target(target_pos)` that returns the power density and phase
                        of the radiation at the given target position.
     - target_pos (tuple or list): The target position where the radiation power density and phase
                                   are to be calculated. The format and dimensionality (e.g., 2D, 3D)
                                   should be consistent with what the `radiation_at_target` method expects.

     Returns:
     - tuple: A tuple containing two lists:
         1. power_densities_at_target (list): A list of power densities at the target position
                                              from each antenna.
         2. phases_at_target (list): A list of phases at the target position from each antenna.

     Note:
     Ensure that the antenna objects in the `antennas` list have the `radiation_at_target` method
     implemented and that it returns the power density and phase in the expected format.
     """
    power_densities_at_target, phases_at_target = [], []
    for antenna in antennas:
        power_density_at_target, phase_at_target = antenna.radiation_at_target(target_pos)
        power_densities_at_target.append(power_density_at_target)
        phases_at_target.append(phase_at_target)

    return power_densities_at_target, phases_at_target

def sum_radiation(power_densities_at_target, phases_at_target):
    """
    Sum up the radiation from multiple antennas at the target position.

    Parameters:
    - powers_at_target: An array of power density at the target by each antenna
    - phases_at_target: An array of phase at the target by each antenna

    Returns:
    - Total power density at the target position.
    """
    E_values = np.sqrt(power_densities_at_target) * np.exp(1j * phases_at_target)
    electric_field = np.sum(E_values)
    total_power_density = np.abs(electric_field) ** 2

    return total_power_density

def antenna_radiation_at_target(antenna, target_pos):
    return antenna.radiation_at_target(target_pos)

class IsotropicAntenna:
    def __init__(self, position, start_phase, pulse_power):
        self.position = position # Position of the antenna
        self.start_phase = start_phase # The phase it start to transmit with
        self.pulse_power = pulse_power  # Power of antenna pulse

    def radiation_at_target(self, target_position):
        """
        Calculate the power density and phase of the wave at the target position.

        Parameters:
        - P: Total power radiated by the isotropic antenna (in watts).
        - r: Distance from the antenna to the target (in meters).
        - lambda_: Wavelength of the radiated wave (in meters).
        - phi_0: Initial phase of the wave at the source (in radians). Default is 0.

        Returns:
        - Power density at the target position (in watts per square meter).
        - Phase of the wave at the target position (in radians).
        """
        r = distance(self.position, target_position)
        power_density_on_target = self.pulse_power / (4 * np.pi * r ** 2)
        phase_on_target = k * r + self.start_phase
        return power_density_on_target, phase_on_target


if __name__ == '__main__':
    target_position = (0, 0, 0)
    pos1 = (1, 0, 0)
    pos2 = (0, 1, 0)
    a1 = IsotropicAntenna(pos1, 0, Pt_antenna)
    a2 = IsotropicAntenna(pos2, 0, Pt_antenna)

    power_densities_at_target, phases_at_target = calculate_radiation_at_target(np.array([a1, a2]), target_position)
    total_power_density = sum_radiation(np.array(power_densities_at_target), np.array(phases_at_target))

    a = 5