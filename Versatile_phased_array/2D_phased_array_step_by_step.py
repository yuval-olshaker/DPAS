from Versatile_phased_array.utils import *


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
    """
    Represents an antenna system.

    This class models an isotropic antenna system with the ability to introduce position and time errors during initialization.

    Parameters:
        position (tuple): A tuple representing the position (x, y, z) of the antenna.
        start_phase (float): The phase at which the antenna starts transmitting.
        pulse_power (float): The power of the antenna pulse.
        pos_error (bool, optional): If True, random position errors are added during initialization.
            Default is True.
        time_error (bool, optional): If True, clock error is simulated by adding a phase error based on
            a normal distribution. Default is True.

    Attributes:
        position (tuple): The position (x, y, z) of the antenna.
        start_phase (float): The phase at which the antenna starts transmitting.
        pulse_power (float): The power of the antenna pulse.
    """

    def __init__(self, position, start_phase, pulse_power, pos_error=True, time_error=True):
        self.desired_position = position  # Position of the antenna
        self.desired_start_phase = start_phase  # The phase it starts to transmit with
        self.pulse_power = pulse_power  # Power of antenna pulse
        if pos_error:  # If we want position error we add it here
            self.position = tuple(position + np.random.normal(0, pos_error_sigma_per_axis, 3))
        else:
            self.position = position
        if time_error:  # If we want clock error we add it here
            # we act like the time error affects the pase as if we changed the range by the error times speed of light
            self.start_phase = start_phase + np.random.normal(0, time_error_sigma) * c * k
        else:
            self.start_phase = start_phase


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
        phase_on_target = normalize_angle(k * target_range + self.start_phase)
        return power_density_on_target, phase_on_target

    def __eq__(self, other):
        if not isinstance(other, IsotropicAntenna):
            return NotImplemented
        return (self.position == other.position and
                self.start_phase == other.start_phase and
                self.pulse_power == other.pulse_power)


def create_isotropic_antennas_planar_array(center):
    """
    Create a planar array of isotropic antennas with positions determined by the grid spacing `dy` and `dz`.

    The function generates positions for each antenna in a planar grid defined by the spacings `dy` and `dz`.
    The entire planar array is then centered around the provided `center` point. Each antenna is an instance
    of the `IsotropicAntenna` class. The x-coordinate for each antenna is fixed at 0.

    Parameters:
    - center (tuple): A tuple of (x, y, z) coordinates representing the center point where the planar array should be positioned.

    Returns:
    - list: A list of `IsotropicAntenna` objects positioned in a planar grid centered around the provided `center` point.
    """
    planar_center = (0, Ny * dy / 2, Nz * dz / 2)
    antennas_array, positions = [], []
    for i in range(Ny + 1):
        for j in range(Nz + 1):
            # we subtract the array center so it will be centered
            x_pos = 0 - planar_center[0] + center[0]
            y_pos = i * dy - planar_center[1] + center[1]
            z_pos = j * dz - planar_center[2] + center[2]
            a_pos = (x_pos, y_pos, z_pos)
            positions.append(a_pos)
            antennas_array.append(IsotropicAntenna(a_pos, 0, Pt_antenna))
    return np.array(antennas_array), positions


def create_random_isotropic_antennas_array(center):
    """
    Create an array of isotropic antennas with random positions within a cube.

    The function generates random positions for each antenna within a cube of side length `array_side_size`.
    The cube is then centered around the provided `center` point. Each antenna is an instance of the
    `IsotropicAntenna` class.

    Parameters:
    - center (tuple): A tuple of (x, y, z) coordinates representing the center point where the cube should be positioned.

    Returns:
    - list: A list of `IsotropicAntenna` objects with random positions within the cube.
    """
    antennas_array = []
    positions = generate_random_3D_points((Ny + 1) * (Nz + 1), center, radius_for_random_antenna_positions)
    for pos in positions:
        antennas_array.append(IsotropicAntenna(pos, 0, Pt_antenna))
    return np.array(antennas_array), positions


def check_array_on_target(antennas_array, target_position):
    """
    Calculate the total power density at a target position from an array of given isotropic antennas.

    The function calculates the power densities and phases at the target position
    for each antenna. Finally, it computes the total power density at the target position by summing
    the contributions from all antennas.

    Parameters:
    - antennas_array: The array of antennas
    - target_position (tuple or list): The target position where the total power density is to be calculated.
                                       The format and dimensionality (e.g., 2D, 3D) should be consistent
                                       with what the `calculate_radiation_at_target` and `IsotropicAntenna`
                                       class expect.

    Returns:
    - float: The total power density at the target position.
    """
    power_densities_at_target, phases_at_target = calculate_radiation_at_target(antennas_array, target_position)
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

    # Calculate the positions using vectorized operations
    x = point[0] + dis * np.cos(ele_mesh) * np.cos(azi_mesh)
    y = point[1] + dis * np.cos(ele_mesh) * np.sin(azi_mesh)
    z = point[2] + dis * np.sin(ele_mesh)

    # Reshape and stack the results
    positions = np.dstack((x, y, z))

    # Reshape and stack the angles
    angles = np.dstack((azi_mesh, ele_mesh))
    return positions, np.degrees(angles)


def calculate_the_density_of_single_antenna_by_given_positions(antenna_pos, target_range):
    """
    Calculate the radiation density of a single isotropic antenna at a given range.

    This function initializes an isotropic antenna at a specified position and calculates
    its radiation density at a given target range.

    Parameters:
    - antenna_pos (tuple or list): The position (x, y, z) of the isotropic antenna.
    - target_range (float): The range (distance) from the antenna at which the radiation
                            density is to be calculated.

    Returns:
    - float: The radiation density of the isotropic antenna at the specified target range.
    """
    antenna = IsotropicAntenna(antenna_pos, 0, Pt_antenna, False, False)
    return antenna.radiation_at_range(target_range)


def compute_phase_shifts(nx, ny, dx, dy, wavelength, theta, phi):
    """
    Computes phase shifts for a 2D phased array.

    Parameters:
    - nx: Number of elements in the x-direction.
    - ny: Number of elements in the y-direction.
    - dx: Spacing between elements in the x-direction.
    - dy: Spacing between elements in the y-direction.
    - wavelength: Wavelength of the signal.
    - theta: Elevation angle in radians.
    - phi: Azimuth angle in radians.

    Returns:
    - 2D numpy array of phase shifts for each element.
    """

    # Create arrays for the x and y indices of the elements
    nx_indices = np.arange(nx)
    ny_indices = np.arange(ny)

    # Meshgrid to create 2D arrays for the x and y indices
    nxx, nyy = np.meshgrid(nx_indices, ny_indices)

    # Compute the phase shift for each element
    delta_phi = -2 * np.pi * ((nxx * dx * np.sin(theta) * np.cos(phi) +
                               nyy * dy * np.sin(theta) * np.sin(phi)) / wavelength)

    return delta_phi


def compute_phase_shifts_for_given_antenna_positions(antenna_positions, theta, phi):
    """
    Computes phase shifts for antennas with random positions.

    Parameters:
    - antenna_positions: A numpy array of shape (N, 3), where N is the number of antennas, and each row is the (x, y, z) position of an antenna.
    - theta: Elevation angle in radians.
    - phi: Azimuth angle in radians.

    Returns:
    - 1D numpy array of phase shifts for each antenna.

     Note:
    - The function assumes a global variable 'wavelength' representing the wavelength of the signal.
    - The phase shifts are computed based on the principle of phase array beamforming, where adjusting
      the phase of signals from individual antennas can steer the direction of the main beam.
    """

    # Unit direction vector for the desired direction
    direction = np.array([np.sin(theta) * np.cos(phi),
                          np.sin(theta) * np.sin(phi),
                          np.cos(theta)])

    # Compute the dot product between the antenna positions and the direction vector
    distance_differences = np.dot(antenna_positions, direction)

    # Normalize by the wavelength and convert to phase shifts
    delta_phases = -2 * np.pi * distance_differences / wavelength

    return delta_phases


def shift_phases_of_antennas_array(antennas_array, theta, phi):
    """
    Shift the phases of antennas in the array based on their positions and given angles.

    This function calculates the phase shifts for each antenna in the array based on their positions
    and the provided angles (theta and phi). It then updates the 'start_phase' attribute of each antenna
    with the computed phase shift. We work on the desired positions - without the error

    Parameters:
    - antennas_array (numpy array): An array of antenna objects. Each antenna object should have
                                 a 'position' attribute representing its location and a 'start_phase'
                                 attribute representing its initial phase.
    - theta (float): The theta angle (in radians) used for calculating the phase shift.
    - phi (float): The phi angle (in radians) used for calculating the phase shift.

    Returns:
    - numpy array: The updated antennas array with modified 'start_phase' attributes.

    Dependencies:
    - This function relies on the following external definitions:
      1. compute_phase_shifts_for_random_positions: A function that calculates the phase shifts
                                                    for given positions and angles.
    """

    # Extract the 'position' attribute using numpy's vectorized operations
    antenna_desired_positions = np.transpose(np.vectorize(lambda antenna: antenna.desired_position)(antennas_array))
    delta_phases = compute_phase_shifts_for_given_antenna_positions(antenna_desired_positions, theta, phi)
    for antenna, delta_phase in zip(antennas_array, delta_phases):
        antenna.start_phase += normalize_angle(delta_phase)
    return antennas_array


def print_max_gain_and_angles():
    """
    Prints the maximum gain value and its corresponding angles from the global 'gains_in_db' array.

    This function searches for the maximum gain value in the 'gains_in_db' array and identifies its
    corresponding angles from the 'angles' array. It then prints the maximum gain and the associated
    azimuth (Phi) and elevation (Theta) angles.

    Global Variables:
    - gains_in_db (numpy array): A 2D array containing gain values in decibels.
    - angles (numpy array): A 2D array containing tuples of (Phi, Theta) angles corresponding to
                            the gain values in 'gains_in_db'.

    Returns:
    - None: This function does not return any value; it only prints the results.
    """

    # Find the maximum value
    max_gain = np.amax(gains_in_db)
    # Find the linear index of the maximum value
    linear_index = np.argmax(gains_in_db)
    # Convert the linear index to 2D indices
    row, col = np.unravel_index(linear_index, gains_in_db.shape)
    # Get the angles
    Phi, Theta = angles[row, col]

    # prints
    print(f"Maximum gain: {max_gain}")
    print(f"Angles of maximum value: (Phi {Phi}, Theta {Theta})")
    return max_gain


if __name__ == '__main__':
    """
    This script simulates the radiation pattern of a planar antenna array with isotropic antennas.

    Key Variables:
    - array_center_pos: The center position of the antenna array.
    - antennas_array: The array of antennas.
    - target_range: The range to the target.
    - positions: The positions of the antennas.
    - angles: The angles corresponding to the positions.
    - density_of_single_antenna: The density of a single antenna at the given position.
    - densities: The densities of the antennas.
    - gains: The gains of the antennas.
    - gains_in_db: The gains of the antennas in decibels.

    Functions:
    - create_antennas_array(): Creates an array of antennas.
    - shift_phases_of_antennas_array(): Shifts the phases of the antennas in the array.
    - generate_positions_at_distance_angles_from_point(): Generates positions at a specified distance and angles from a point.
    - calculate_the_density_of_single_antenna_by_given_positions(): Calculates the density of a single antenna at the given positions.
    - check_array_on_target(): Checks the array on the target and returns the density.
    - print_max_gain_and_angles(): Prints the maximum gain and the corresponding angles.

    The script:
    1. Initializes the center position of the antenna array.
    2. Creates an array of antennas and shifts their phases.
    3. Generates positions and angles for the antennas at a specified distance from the center position.
    4. Calculates the density of a single antenna at the given positions.
    5. Iterates over the positions and calculates the density and gain for each position.
    6. Converts the gains to decibels.
    7. Plots the radiation pattern of the antenna array.
    """
    array_center_pos = (0, 0, 0)
    antennas_array, antennas_positions = create_isotropic_antennas_planar_array(array_center_pos)
    # antennas_array, antennas_positions = create_random_isotropic_antennas_array(array_center_pos)
    antennas_array = shift_phases_of_antennas_array(antennas_array, np.radians(100), np.radians(10))
    # print(centroid_of_convex_hull(antennas_positions))
    target_range = 10000
    positions, angles = generate_positions_at_distance_angles_from_point(array_center_pos, PHI, THETA, target_range)
    density_of_single_antenna, _ = calculate_the_density_of_single_antenna_by_given_positions(array_center_pos, target_range)
    densities = np.zeros((positions.shape[0], positions.shape[1]))
    gains = np.zeros((positions.shape[0], positions.shape[1]))
    for i in range(positions.shape[0]):
        for j in range(positions.shape[1]):
            den = check_array_on_target(antennas_array, positions[i, j])
            gain = den / density_of_single_antenna
            densities[i, j] = den
            gains[i, j] = gain
            # print(f"Density at Position: {positions[i, j]}, Angles (Azimuth, Elevation): {angles[i,j]} is: {den}, the gain is {gain}")

    gains_in_db = 10 * np.log10(gains)
    max_gain = print_max_gain_and_angles()

    # Plot
    PHI_deg = np.rad2deg(PHI)
    THETA_deg = np.rad2deg(THETA)
    plt.figure()
    pcm = plt.imshow(gains_in_db, extent=[PHI_deg.min(), PHI_deg.max(), THETA_deg.min(), THETA_deg.max()],
                     origin='lower', aspect='auto', cmap='viridis', vmin=-max_gain, vmax=max_gain)
    plt.colorbar(label='Radiation Pattern (dB)')
    plt.title('Planar Antenna Array with Isotropic Antennas Radiation Pattern')
    plt.xlabel('Phi')
    plt.ylabel('Theta')

    # Add tooltips with gain values using mplcursors
    cursor = mplcursors.cursor(hover=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(f"Gain (dB): {gains_in_db[sel.target.index]:.2f}"))

    plt.show()

    # for pos in positions:
    #     den = check_array_on_target(pos)
    #     densities.append(den)
    #     print('density at pos: ' + str(pos) + ' is: ' + str(den))

    a = 5
