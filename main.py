from rocket import Rocket, Motor, Fins, NoseCone
from flight import Flight

"""
Main function for running the simulation
"""
def main() -> None:
    # Rocket components
    motor = Motor(
        mass=4.766,  # kg
        length=0.795,  # m
        peak_thrust=3415,  # Newtons
        burn_time=2.9  # Seconds
    )
    nose_cone = NoseCone(
        mass=1.721,  # kg
        diameter=0.156,  # m
        length=0.776  # m
    )
    fins = Fins(
        position=0.01,  # m
        fin_mass=1.369,  # kg
        num_fins=4,
        width=0.003,  # m
        root_chord=0.3,  # m
        tip_chord=0.028,  # m
        semi_span=0.14,  # m
        sweep_length=0.225  # m
    )

    # Rocket configuration
    rocket = Rocket(
        motor=motor,
        body_mass=15.238,  # kg
        diameter=0.156,  # m
        body_length=1.5  # m
    )
    rocket.attach_nose_cone(nose_cone)
    rocket.attach_fins(fins)

    # Run the simulation (without drag data)
    flight = Flight(rocket)
    flight.summary()
    flight.write_to_file("Data/sim_data.csv")
    flight.plot("times", ("thrust_forces", "Thrust"), ("drag_forces", "Drag"), ("gravitational_forces", "Weight"), title="Forces Plot", x_label="Time (s)", y_label="Force (N)")
    flight.plot("times", ("altitudes", "Altitude"), ("velocities", "Velocity (m/s)"), ("accelerations", "Acceleration (m/s^2)"), title="Altitude Plot", x_label="Time (s)", y_label="Altitude")

    # Run the simulation (with drag data)
    flight = Flight(rocket, drag_data="Data/drag_data.csv")
    flight.summary()
    flight.write_to_file("Data/sim_data_cfd.csv")
    flight.plot("times", ("thrust_forces", "Thrust"), ("drag_forces", "Drag"), ("gravitational_forces", "Weight"), title="Forces Plot", x_label="Time (s)", y_label="Force (N)")
    flight.plot("times", ("altitudes", "Altitude"), ("velocities", "Velocity (m/s)"), ("accelerations", "Acceleration (m/s^2)"), title="Altitude Plot", x_label="Time (s)", y_label="Altitude")

if __name__ == "__main__":
    main()
