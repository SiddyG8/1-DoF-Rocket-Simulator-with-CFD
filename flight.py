from rocket import Rocket
import constants as c
from numpy import e, sin, log10

class Flight:
    def __init__(self, rocket: Rocket, dt: float) -> None:
        self.rocket = rocket
        self.dt = dt
        self.t = 0
        self.z = 0
        self.vz = 0
        self.az = 0
        self.phase = "pad"

    @property
    def altitude(self) -> float:
        return self.z

    def drag_coefficient(self) -> float:
        return pow(e, -1.2 * self.M) * sin(self.M) + (self.M / 6) * log10(self.M + 1)

    def mach_number(self) -> float:
        pass

    def simulate(self) -> None:
        pass

if __name__ == "__main__":
    pass
