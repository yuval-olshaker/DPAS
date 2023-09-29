from phased_arrays.utils import *

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
        watts_on_target = self.pulse_power / (4 * np.pi * r ** 2)
        phase_on_target = k * r + self.start_phase
        return watts_on_target, phase_on_target


if __name__ == '__main__':
    a = 5