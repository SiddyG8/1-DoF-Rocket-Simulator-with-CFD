from motor import Motor

class Rocket:
    def __init__(self, motor, wet_mass, fuel_mass) -> None:
        self.motor = motor
        self.wet_mass = wet_mass
        self.fuel_mass = fuel_mass

    @property
    def dry_mass(self) -> float:
        return self.wet_mass - self.fuel_mass

if __name__ == "__main__":
    pass
