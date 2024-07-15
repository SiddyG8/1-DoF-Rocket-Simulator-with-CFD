from rocket import Rocket, Motor, Fins, NoseCone, Parachute
from flight import Flight


# Project TODOs:
# Add docstrings
# Update README.md

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
    main_parachute = Parachute(
        diameter=2.44,  # m
        drag_coefficient=1.655,
        mass=1.2  # kg
    )
    drogue_parachute = Parachute(
        diameter=0.914,  # m
        drag_coefficient=0.655,
        mass=0.8  # kg
    )

    # Rocket configuration
    rocket = Rocket(
        motor=motor,
        body_mass=13.238,  # kg
        diameter=0.156,  # m
        body_length=1.5  # m
    )
    rocket.attach_nose_cone(nose_cone)
    rocket.attach_fins(fins)
    rocket.attach_main_parachute(main_parachute)
    rocket.attach_drogue_parachute(drogue_parachute)

    # Run the simulation (without drag data)
    flight = Flight(rocket, deploy_altitude=457)
    flight.summary()
    flight.write_to_file("Data/sim_data.csv")
    flight.plot("times", ("thrust_forces", "Thrust"), ("drag_forces", "Drag"), ("gravitational_forces", "Weight"), title="Forces Plot", x_label="Time (s)", y_label="Force (N)")
    flight.plot("times", ("altitudes", "Altitude"), ("velocities", "Velocity (m/s)"), ("accelerations", "Acceleration (m/s^2)"), title="Altitude Plot", x_label="Time (s)", y_label="Altitude")

    # Run the simulation (with drag data)
    flight = Flight(rocket, deploy_altitude=457, drag_data="Data/drag_data.csv")
    flight.summary()
    flight.write_to_file("Data/sim_data_cfd.csv")
    flight.plot("times", ("thrust_forces", "Thrust"), ("drag_forces", "Drag"), ("gravitational_forces", "Weight"), title="Forces Plot (CFD)", x_label="Time (s)", y_label="Force (N)")
    flight.plot("times", ("altitudes", "Altitude"), ("velocities", "Velocity (m/s)"), ("accelerations", "Acceleration (m/s^2)"), title="Altitude Plot (CFD)", x_label="Time (s)", y_label="Altitude")

if __name__ == "__main__":
    main()
