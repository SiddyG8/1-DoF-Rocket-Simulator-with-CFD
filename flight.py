from rocket import Rocket

class Flight:
    def __init__(self, rocket, dt) -> None:
        self.rocket = rocket
        self.dt = dt
        self.t = 0
        self.z = 0
        self.vz = 0
        self.az = 0
        self.phase = "On pad"

    def simulate(self) -> None:
        """
        Simulate the flight of the rocket.
        """
        pass

if __name__ == "__main__":
    pass
