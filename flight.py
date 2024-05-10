from rocket import Rocket
import constants as c
from numpy import e, sin, log10

class Flight:
    def __init__(self, rocket: Rocket, dt: float = 0.05) -> None:
        self.rocket = rocket
        self.dt = dt
        self.t = 0
        self.z = []
        self.vz = []
        self.az = []
        self.phase = []

    @property
    def altitude(self) -> list[float]:
        return self.z

    def temperature(self) -> float:
        return c.T0 - c.L * self.z

    def pressure(self) -> float:
        return c.p0 * (self.temperature / c.T0) ** (-c.g0 / (c.R * c.L))

    def drag_coefficient(self) -> float:
        return pow(e, -1.2 * self.M) * sin(self.M) + (self.M / 6) * log10(self.M + 1)

    def drag_force(self) -> float:
        return 0.5 * c.rho * self.vz * self.vz * self.rocket.wetted_area * self.drag_coefficient()

    def simulate(self) -> None:
        pass

if __name__ == "__main__":
    pass
