import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import mplcursors
import plotly.graph_objects as go

# Parameters
frequency = 1e9  # 1 GHz
c = 3e8  # Speed of light
wavelength = c / frequency # 0.3 m - 30 cm
k = 2 * np.pi / wavelength #wave number
epsilon = 1e-10
antenna_size = wavelength / 2


# Array dimensions
Ny = 10  # Number of elements along y
Nz = Ny  # Number of elements along z
dy = 0.5 * wavelength  # Spacing along y
dz = 0.5 * wavelength  # Spacing along z
x = np.arange(Ny) * dy
y = np.arange(Nz) * dz
radius_for_random_antenna_positions = 10 # 10m radius

# Angles
azimuth_range = np.radians(30) # 120 degrees
elevation_range = np.radians(30) # 60 degrees
PHI = np.linspace(- azimuth_range / 2, azimuth_range / 2, 50)
THETA = np.linspace(-elevation_range / 2, elevation_range / 2, 50)
azimuth_resolution = 1 # milliradian
azimuth_resolution = (azimuth_resolution / (1000 * 2 * np.pi)) * 360 # degrees
elevation_resolution = 2 # milliradian
elevation_resolution = (elevation_resolution / (1000 * 2 * np.pi)) * 360 # degrees

# Pulse energy at each antenna
Pt_antenna = 1 # in watt per pulse
Pt_array = Pt_antenna * Ny * Nz

# Target
R = 100 * 1e3
RCS = 10

# Noise
temperature = 290
kB = 1.38 * 1e-23
pulse_bandwidth = 50 * 1e6

# Params for SNR
G = 5000 # 37 DB
Ae = Ny * Nz * 1 # 1 meter spacing for antenna size. affective size is the whole size

# Errors
pos_error_sigma_per_axis = 0 / 100 # 0.5 cm error per axis

def distance(pos1, pos2):
    """
        Calculate the Euclidean distance between two points.

        Parameters:
        - pos1 (tuple or list): Coordinates of the first point. Can be in N-dimensional space.
        - pos2 (tuple or list): Coordinates of the second point. Must have the same dimensionality as pos1.

        Returns:
        - float: The Euclidean distance between pos1 and pos2.
    """
    return np.linalg.norm(np.array(pos1) - np.array(pos2))

def normalize_angle(angle):
    """
    Normalize an angle to the range [-pi, pi].

    Parameters:
    - angle (float): The angle in radians.

    Returns:
    - float: The normalized angle in the range [-pi, pi].
    """
    normalized_angle = angle % (2 * np.pi)
    if normalized_angle > np.pi:
        normalized_angle -= 2 * np.pi
    return normalized_angle

def centroid_of_convex_hull(points):
    """
    Calculate the centroid of the convex hull formed by an array of 3D points.

    Parameters:
    - points: Array of 3D points.

    Returns:
    - Centroid of the convex hull.
    """
    hull = ConvexHull(points)
    total_volume = 0
    centroid_sum = np.zeros(3)

    for simplex in hull.simplices:
        # Form a tetrahedron using the reference point and the face
        tetra = [hull.points[simplex[0]], hull.points[simplex[1]], hull.points[simplex[2]], hull.points[0]]
        volume = ConvexHull(tetra).volume
        total_volume += volume

        # Calculate the centroid of the tetrahedron
        tetrahedron_centroid = np.mean(tetra, axis=0)
        centroid_sum += tetrahedron_centroid * volume

    # Calculate the final centroid
    centroid = centroid_sum / total_volume
    return centroid


def centroid_of_planar_convex_hull(points):
    """
    Calculate the centroid of the convex hull for points that all lie on the same plane.

    Parameters:
    - points: Array of 3D points.

    Returns:
    - Centroid of the convex hull.
    """
    points = np.array(points)
    # Project points to 2D by removing the constant coordinate
    points_2d = points[:, 1:]

    # Compute the convex hull of the 2D points
    hull = ConvexHull(points_2d)

    # Calculate the centroid of the 2D convex hull
    centroid_2d = np.mean(hull.points[hull.vertices], axis=0)

    # Map the 2D centroid back to 3D
    centroid = np.array([points[0, 0], centroid_2d[0], centroid_2d[1]])

    return centroid


def generate_random_3D_points(n, centroid, r):
    """
    Generate n random 3D points around a specified centroid within a radius r.

    Parameters:
    - n: Number of random points to generate.
    - centroid: Tuple (x_centroid, y_centroid, z_centroid).
    - r: Radius within which to generate random points.

    Returns:
    - Array of shape (n, 3) containing the random 3D points.
    """
    points = []
    for _ in range(n):
        # Generate random direction (unit vector)
        theta = 2 * np.pi * np.random.rand()
        phi = np.arccos(2 * np.random.rand() - 1)
        x_dir = np.sin(phi) * np.cos(theta)
        y_dir = np.sin(phi) * np.sin(theta)
        z_dir = np.cos(phi)

        # Scale by random length between 0 and r
        length = r * np.cbrt(np.random.rand())  # Cube root for uniform distribution in 3D
        x_offset = x_dir * length
        y_offset = y_dir * length
        z_offset = z_dir * length

        # Add offset to centroid
        x_point = centroid[0] + x_offset
        y_point = centroid[1] + y_offset
        z_point = centroid[2] + z_offset

        points.append([x_point, y_point, z_point])

    return np.array(points)
