from rocket import Rocket
from flight import Flight
import matplotlib.pyplot as plt


def main() -> None:
    rocket = Rocket(
        mass=27.201, # kg
        fuel_mass=4.766, # kg
        diameter=0.156, # meters
        T0=3415, # Newtons
        tb=2.9 # Seconds
    )
    flight = Flight(rocket)
    # print(max(flight.z))
    # plt.plot(flight.t, flight.vz)
    plt.plot(flight.t, flight.thrust_forces, label="Thrust")
    plt.plot(flight.t, flight.drag_forces, label="Drag")
    plt.plot(flight.t, flight.gravitional_forces, label="Gravitation")
    print(max(flight.z))
    print(max(flight.drag_forces))
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
