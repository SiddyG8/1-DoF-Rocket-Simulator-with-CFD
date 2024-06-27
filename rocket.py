from dataclasses import dataclass
from motor import Motor
from nose_cone import NoseCone
from fins import Fins
from math import pi

@dataclass
class Rocket:
    mass: float
    fuel_mass: float
    diameter: float
    motor: Motor
    nose_cone: NoseCone
    fins: Fins

    @property
    def delta_mass(self) -> float:
        return -self.fuel_mass / self.motor.tb

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2

    def thrust(self, t) -> float:
        return self.motor.thrust(t)


if __name__ == "__main__":
    pass
