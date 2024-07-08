from dataclasses import dataclass
from math import pi, e, log


@dataclass
class Motor:
    mass: float
    length: float
    peak_thrust: float
    burn_time: float

    @property
    def delta_mass(self) -> float:
        return self.mass / self.burn_time

    def thrust(self, t) -> float:
        # Equation driven thrust curve that estimates a motor's thrust over time based on its peak thrust and burn time
        return max(0, self.peak_thrust * (1 - (10**-5) * (e**(log(10**5)/self.burn_time * t))))


"""
Represents the (Von Kármán) nose cone of a rocket.
"""
class NoseCone:
    X_N = 0.466  # Nose cone length ratio (Ogive)

    def __init__(self, mass: float, diameter: float, length: float) -> None:
        self.mass = mass
        self.diameter = diameter
        self.length = length

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2


@dataclass
class Fins:
    position: float
    fin_mass: float
    num_fins: int
    width: float
    root_chord: float
    tip_chord: float
    semi_span: float
    sweep_length: float

    @property
    def total_mass(self) -> float:
        return self.num_fins * self.fin_mass

    @property
    def wetted_area(self) -> float:
        return self.num_fins * self.width * self.semi_span


class Rocket:
    def __init__(self, motor: Motor, body_mass: float, diameter: float, body_length: float) -> None:
        self.body_mass = body_mass
        self.diameter = diameter
        self.body_length = body_length
        self.motor = motor
        self.nose_cone = None
        self.fins = None

    @property
    def total_mass(self) -> float:
        total_mass = self.body_mass

        if self.motor:
            total_mass += self.motor.mass
        if self.nose_cone:
            total_mass += self.nose_cone.mass
        if self.fins:
            total_mass += self.fins.total_mass

        return total_mass

    @property
    def wetted_area(self) -> float:
        wetted_area = 0

        if self.nose_cone:
            wetted_area += self.nose_cone.wetted_area
        if self.fins:
            wetted_area += self.fins.wetted_area

        return wetted_area

    def attach_nose_cone(self, nose_cone: NoseCone) -> None:
        self.nose_cone = nose_cone

    def attach_fins(self, fins: Fins) -> None:
        self.fins = fins

    def thrust(self, t) -> float:
        return self.motor.thrust(t)


if __name__ == "__main__":
    pass
