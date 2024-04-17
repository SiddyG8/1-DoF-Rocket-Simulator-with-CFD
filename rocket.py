from motor import Motor

class Rocket:
    def __init__(self, motor, wet_mass, fuel_mass, diameter, fins = None) -> None:
        self.motor = motor
        self.fins = fins
        self.wet_mass = wet_mass
        self.fuel_mass = fuel_mass
        self.diameter = diameter

    @property
    def dry_mass(self) -> float:
        return self.wet_mass - self.fuel_mass

if __name__ == "__main__":
    pass
