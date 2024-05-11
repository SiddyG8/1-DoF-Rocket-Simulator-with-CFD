import math

class Rocket:
    def __init__(self, motor, wet_mass, fuel_mass, diameter, fins = None, nose_cone = None) -> None:
        self.motor = motor
        self.fins = fins
        self.nose_cone = nose_cone
        self.wet_mass = wet_mass
        self.fuel_mass = fuel_mass
        self.diameter = diameter

    @property
    def dry_mass(self) -> float:
        return self.wet_mass - self.fuel_mass
    
    @property
    def wetted_area(self) -> float:
        wetted_area = math.pi * (self.diameter / 2) ** 2
        if self.nose_cone and self.nose_cone.diameter > self.diameter:
            wetted_area = self.nose_cone.wetted_area
        if self.fins:
            wetted_area += self.fins.wetted_area
        return wetted_area


if __name__ == "__main__":
    pass
