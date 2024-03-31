from rocket import Rocket

class Flight:
    def __init__(self, rocket, dt) -> None:
        self.rocket = rocket
        self.dt = dt
        self.t = 0

    def simulate(self) -> None:
        """
        Simulate the flight of the rocket.
        """
        self.rocket.motor.burnout()
        while self.t < self.rocket.motor.tb:
            self.rocket.motor.thrust()
            self.rocket.accelerate()
            self.rocket.propagate()
            self.t += self.dt

if __name__ == "__main__":
    pass
