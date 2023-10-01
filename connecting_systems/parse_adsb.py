"""
    Chat GPT generated code that pars a massage from the ICD PDF of ADSB protocol.
    Used a plugin that can parse PDFs.
    The parsing is not perfect for other massages but the ability exists.
"""

class AircraftOperationalStatus:
    def __init__(self, binary_data):
        self.data = binary_data

    def get_TCAS_resolution_advisory(self):
        return "TCAS RA active" if self.data[0] == '1' else "TCAS II or ACAS RA not active"

    def get_target_trajectory_change_capability(self):
        mapping = {
            '00': "no capability for Trajectory Change Reports",
            '01': "support for TC+0 reports only",
            '10': "support for multiple TC reports",
            '11': "reserved"
        }
        return mapping[self.data[1:3]]

    def get_target_state_report_capability(self):
        return "capable of supporting target State Reports" if self.data[4] == '1' else "no capability to support Target State Reports"

    def get_air_referenced_velocity_capability(self):
        return "capable of generating ARV-reports" if self.data[5] == '1' else "no capability to generate ARV-reports"

    def get_CDTI_operational_status(self):
        return "CDTI operational" if self.data[6] == '1' else "CDTI not operational"

    def get_TCAS_system_status(self):
        return "TCAS operational" if self.data[7] == '0' else "TCAS not operational"

    def get_antenna_status(self):
        return "Single Antenna only" if self.data[8] == '1' else "Antenna Diversity"

    def display_status(self):
        print(f"TCAS Resolution Advisory: {self.get_TCAS_resolution_advisory()}")
        print(f"Target Trajectory Change Capability: {self.get_target_trajectory_change_capability()}")
        print(f"Target State Report Capability: {self.get_target_state_report_capability()}")
        print(f"Air-Referenced Velocity Capability: {self.get_air_referenced_velocity_capability()}")
        print(f"CDTI Operational Status: {self.get_CDTI_operational_status()}")
        print(f"TCAS System Status: {self.get_TCAS_system_status()}")
        print(f"Antenna Status: {self.get_antenna_status()}")


def main():
    # Simulated binary data for Aircraft Operational Status
    binary_data = "010101010"
    status = AircraftOperationalStatus(binary_data)
    status.display_status()


if __name__ == "__main__":
    main()
