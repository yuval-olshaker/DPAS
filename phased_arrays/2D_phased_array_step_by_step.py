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
    """
    Calculates the radiation_at_target of the antenna and target given
    :param antenna: the antenna
    :param target_pos: the position of the target
    :return: antenna.radiation_at_target(target_pos)
    """
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
        - target_position: Position of the target.

        Returns:
        - Power density at the target position (in watts per square meter).
        - Phase of the wave at the target position (in radians).
        """
        target_range = distance(self.position, target_position)
        return self.radiation_at_range(target_range)

    def radiation_at_range(self, target_range):
        """
        Calculate the power density and phase of the wave at the target position.

        Parameters:
        - P: Total power radiated by the isotropic antenna (in watts).
        - target_range: range to the target.

        Returns:
        - Power density at the target position (in watts per square meter).
        - Phase of the wave at the target position (in radians).
        """
        power_density_on_target = self.pulse_power / (4 * np.pi * target_range ** 2)
        phase_on_target = k * target_range + self.start_phase
        return power_density_on_target, phase_on_target


def check_array_on_target(target_position):
    """
        Calculate the total power density at a target position from an array of isotropic antennas.

        The function initializes an array of isotropic antennas based on predefined grid dimensions (Nx, Ny)
        and spacing (dx, dy). It then calculates the power densities and phases at the target position
        for each antenna. Finally, it computes the total power density at the target position by summing
        the contributions from all antennas.

        Parameters:
        - target_position (tuple or list): The target position where the total power density is to be calculated.
                                           The format and dimensionality (e.g., 2D, 3D) should be consistent
                                           with what the `calculate_radiation_at_target` and `IsotropicAntenna`
                                           class expect.

        Returns:
        - float: The total power density at the target position.

        Dependencies:
        - This function relies on the following external definitions:
            1. Nx, Ny: Grid dimensions for the antenna array.
            2. dx, dy: Spacing between antennas in the x and y directions.
            3. Pt_antenna: Power of each isotropic antenna.
            4. IsotropicAntenna: A class representing an isotropic antenna. It should have a method or property
                                 to return the radiation at a target position.
            5. calculate_radiation_at_target: A function that calculates power densities and phases at a target
                                              position for a given set of antennas.
            6. sum_radiation: A function that sums up the power densities and phases to compute the total power
                              density at the target position.

        Note:
        Ensure that the dependencies (Nx, Ny, dx, dy, Pt_antenna, IsotropicAntenna, calculate_radiation_at_target,
        and sum_radiation) are properly defined and available in the scope where this function is used.
        """
    antennas = []
    for i in range(Nx + 1):
        for j in range(Ny + 1):
            pos = (0, i * dx, j * dy)
            antennas.append(IsotropicAntenna(pos, 0, Pt_antenna))

    power_densities_at_target, phases_at_target = calculate_radiation_at_target(np.array(antennas), target_position)
    total_power_density = sum_radiation(np.array(power_densities_at_target), np.array(phases_at_target))
    return total_power_density

def generate_positions_at_distance_angles_from_point(point, angles_azi, angles_ele, dis):
    """
    Generate positions at a given distance for all combinations of azimuth and elevation angles from a reference point.

    Parameters:
    - point (tuple or list): The reference point (x, y, z).
    - angles_azi (list or numpy array): List of azimuth angles in radians.
    - angles_ele (list or numpy array): List of elevation angles in radians.
    - dis (float): Distance from the reference point.

    Returns:
    - numpy array: Array of positions (x, y, z) at the given distance for all combinations of azimuth and elevation angles.
    """
    # Create a meshgrid for all combinations
    azi_mesh, ele_mesh = np.meshgrid(angles_azi, angles_ele)

    # Calculate the new positions using vectorized operations
    x = point[0] + dis * np.cos(ele_mesh) * np.cos(azi_mesh)
    y = point[1] + dis * np.cos(ele_mesh) * np.sin(azi_mesh)
    z = point[2] + dis * np.sin(ele_mesh)

    # Reshape and stack the results
    positions = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

    # Reshape and stack the angles
    angles = np.vstack((azi_mesh.ravel(), ele_mesh.ravel())).T

    return positions, np.degrees(angles)

def calculate_the_density_of_single_antenna_by_given_positions(antenna_pos, target_range):
    antenna = IsotropicAntenna(antenna_pos, 0, Pt_antenna)
    return antenna.radiation_at_range(target_range)

if __name__ == '__main__':
    array_center_pos = (0, Nx * dx / 2, Ny * dy / 2)
    target_range = 10000
    positions, angles = generate_positions_at_distance_angles_from_point(array_center_pos, PHI, THETA, target_range)
    density_of_single_antenna, _ = calculate_the_density_of_single_antenna_by_given_positions(array_center_pos, target_range)
    densities, gains = [], []
    for pos, ang in zip(positions, angles):
        den = check_array_on_target(pos)
        gain = den / density_of_single_antenna
        densities.append(den)
        gains.append(gain)
        print(f"Density at Position: {pos}, Angles (Azimuth, Elevation): {ang} is: {den}, the gain is {gain}")

    # for pos in positions:
    #     den = check_array_on_target(pos)
    #     densities.append(den)
    #     print('density at pos: ' + str(pos) + ' is: ' + str(den))

    a = 5