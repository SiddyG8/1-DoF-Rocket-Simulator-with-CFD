from rocket import Rocket
from flight import Flight

class Main:
    def __init__(self) -> None:
        self.rocket = Rocket(
            mass=27.201, # kg
            fuel_mass=4.766, # kg
            diameter=0.156, # meters
            T0=3983, # Newtons
            tb=2.9 # Seconds
        )
        self.flight = Flight(self.rocket)
        self.flight.plot("time", "altitude", "Time (s)", "Altitude (m)")


if __name__ == "__main__":
    Main()
