import math
import matplotlib.pyplot as plt


def radar_horizon(radar_height, target_height):
    # Radius of the Earth in kilometers
    R = 6371

    # Calculate the distance to the horizon for both the radar and the target
    d_radar = math.sqrt(2 * radar_height * R)
    d_target = math.sqrt(2 * target_height * R)

    # Return the sum of both distances
    return d_radar + d_target


# Convert distance to time (in seconds)
def distance_to_time(distance):
    speed = 0.85  # km/s
    return distance / speed

if __name__ == '__main__':

    # Fixed target height
    target_height = 0.015  # 15m in km

    # Varying radar heights from 0 to 5000m
    radar_heights = [i / 1000 for i in range(0, 501)] # in km
    distances = [radar_horizon(radar_height, target_height) for radar_height in radar_heights]
    times = [distance_to_time(distance) for distance in distances]

    # Highlighted radar heights
    highlighted_heights = [0.02, 0.1, 0.15, 0.2, 0.4, 0.5]
    highlighted_distances = [radar_horizon(height, target_height) for height in highlighted_heights]
    highlighted_times = [distance_to_time(distance) for distance in highlighted_distances]

    # Distance graph
    plt.figure(figsize=(10, 6))
    plt.plot(radar_heights, distances, label='Radar Target Visibility', color='blue')
    plt.scatter(highlighted_heights, highlighted_distances, color='red', zorder=5)
    for height, distance in zip(highlighted_heights, highlighted_distances):
        plt.annotate(f"{height * 1000}m", (height, distance), textcoords="offset points", xytext=(0, 10), ha='center')
    plt.xlabel('Radar Height (km)')
    plt.ylabel('Distance (km)')
    plt.title('Radar Target Visibility vs Radar Height - Target height 15m, Speed 850m/s')
    plt.legend()
    plt.grid(True)

    # Add distances of highlighted points to y-axis for distance graph
    current_yticks_distance = plt.yticks()[0]
    updated_yticks_distance = sorted(list(set(current_yticks_distance).union(set(highlighted_distances))))
    plt.yticks(updated_yticks_distance,
               [f"{tick:.2f} km" if tick in highlighted_distances else str(tick) for tick in updated_yticks_distance])

    plt.tight_layout()
    plt.show()

    # Time graph
    plt.figure(figsize=(10, 6))
    plt.plot(radar_heights, times, label='Time from Detection to Impact', color='blue')
    plt.scatter(highlighted_heights, highlighted_times, color='red', zorder=5)
    for height, time in zip(highlighted_heights, highlighted_times):
        plt.annotate(f"{height * 1000}m", (height, time), textcoords="offset points", xytext=(0, 10), ha='center')
    plt.xlabel('Radar Height (km)')
    plt.ylabel('Time (seconds)')
    plt.title('Time from Detection to Impact vs Radar Height - Target height 15m, Speed 850m/s')
    plt.legend()
    plt.grid(True)

    # Add times of highlighted points to y-axis for time graph
    current_yticks_time = plt.yticks()[0]
    updated_yticks_time = sorted(list(set(current_yticks_time).union(set(highlighted_times))))
    plt.yticks(updated_yticks_time,
               [f"{tick:.2f} s" if tick in highlighted_times else str(tick) for tick in updated_yticks_time])

    plt.tight_layout()
    plt.show()
