from math import pi, e, log

class Rocket:
    def __init__(self, mass, fuel_mass, diameter, T0, tb) -> None:
        self.mass = mass
        self.fuel_mass = fuel_mass
        self.diameter = diameter
        self.T0 = T0
        self.tb = tb

    @property
    def delta_mass(self) -> float:
        return -self.fuel_mass / self.tb

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2

    def thrust(self, t) -> float:
        return max(0, self.T0 * (1 - (10**-5) * (e**(log(10**5)/self.tb * t))))

if __name__ == "__main__":
    pass
