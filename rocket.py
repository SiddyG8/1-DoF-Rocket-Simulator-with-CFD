from math import pi, e, log

class Rocket:
    def __init__(self, mass, fuel_mass, diameter, motor) -> None:
        self.mass = mass
        self.fuel_mass = fuel_mass
        self.diameter = diameter
        self.motor = motor

    @property
    def delta_mass(self) -> float:
        return -self.fuel_mass / self.tb

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2

if __name__ == "__main__":
    pass
