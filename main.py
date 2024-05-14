from rocket import Rocket
from flight import Flight


def main() -> None:
    rocket = Rocket(
        mass=27.201, # kg
        fuel_mass=4.766, # kg
        diameter=0.156, # meters
        T0=3415, # Newtons
        tb=2.9 # Seconds
    )
    flight = Flight(rocket)

    # Printing flight statistics
    print(f"Apogee: {max(flight.z):.2f} meters at {flight.t[flight.z.index(max(flight.z))]:.2f} seconds")
    print(f"Max velocity: {max(flight.vz):.2f} m/s at {flight.t[flight.vz.index(max(flight.vz))]:.2f} seconds")
    print(f"Flight time: {flight.t[-1]:.2f} seconds")

    # # Forces plot
    # flight.plot(
    #     "t",
    #     ("thrust_forces", "Thrust"),
    #     ("drag_forces", "Drag"),
    #     ("gravitional_forces", "Gravitation"),
    #     x_label="Time (s)",
    #     y_label="Force (N)"
    # )

    # Altitude plot
    flight.plot(
        "t",
        # ("z", "Altitude"),
        "z",
        # ("vz", "Velocity (m/s)"),
        # "az",
        x_label="Time (s)",
        y_label="Altitude (m)"
    )


if __name__ == "__main__":
    main()
