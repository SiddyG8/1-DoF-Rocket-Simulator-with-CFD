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

    flight = Flight(rocket, drag_data="Data/aggregated_coefficients.csv")

    # Printing flight statistics
    print(f"Apogee: {max(flight.z):.2f} meters at {flight.t[flight.z.index(max(flight.z))]:.2f} seconds")
    print(f"Max velocity: {max(flight.vz):.2f} m/s at {flight.t[flight.vz.index(max(flight.vz))]:.2f} seconds")
    print(f"Flight time: {flight.t[-1]:.2f} seconds")
    print(f"Max Q : {max(flight.dynamic_pressures)/1000:.2f} kPa at {flight.t[flight.dynamic_pressures.index(max(flight.dynamic_pressures))]:.2f} seconds")
    print(f"Starting TWR: {flight.twr[0]:.2f}")
    print(f"Off Rod Velocity: {flight.off_rod_velocity(5)[1]:.2f} m/s at {flight.off_rod_velocity(5)[0]:.2f} seconds")

    # Forces plot
    flight.plot(
        "t",
        ("thrust_forces", "Thrust"),
        ("drag_forces", "Drag"),
        ("gravitational_forces", "Gravitation"),
        x_label="Time (s)",
        y_label="Force (N)"
    )

    # Altitude plot
    flight.plot(
        "t",
        ("z", "Altitude"),
        ("vz", "Velocity (m/s)"),
        ("az", "Acceleration (m/s^2)"),
        x_label="Time (s)",
        y_label="Altitude (m)"
    )

if __name__ == "__main__":
    main()
