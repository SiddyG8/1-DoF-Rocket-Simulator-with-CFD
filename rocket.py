from dataclasses import dataclass
from math import pi, e, log, sqrt


class Motor:
    def __init__(self, position: float, mass: float, peak_thrust: float, burn_time: float) -> None:
        self.position = position
        self.mass = mass
        self.peak_thrust = peak_thrust
        self.burn_time = burn_time
        self._from_curve = False

    @classmethod
    def from_thrust_curve(cls, position, filepath: str) -> None:
        raise NotImplementedError

    @property
    def delta_mass(self) -> float:
        return self.mass / self.burn_time

    def thrust(self, t) -> float:
        if self._from_curve:
            raise NotImplementedError

        # Equation driven thrust curve that estimates a motor's thrust over time based on its peak thrust and burn time
        return max(0, self.peak_thrust * (1 - (10**-5) * (e**(log(10**5)/self.burn_time * t))))


"""
Represents the (Von Kármán) nose cone of a rocket.
"""
class NoseCone:
    X_N = 0.466  # Nose cone length ratio (Ogive)

    def __init__(self, position: float, mass: float, diameter: float, length: float) -> None:
        self.position = position
        self.mass = mass
        self.diameter = diameter
        self.length = length

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2


@dataclass
class Fins:
    position: float
    num_fins: int
    width: float
    root_chord: float
    tip_chord: float
    semi_span: float
    sweep_length: float

    @property
    def wetted_area(self) -> float:
        return self.num_fins * self.width * self.semi_span


class Rocket:
    def __init__(self, mass: float, diameter: float, length: float) -> None:
        self.mass = mass
        self.diameter = diameter
        self.length = length
        self.motor = None
        self.nose_cone = None
        self.fins = None

    @property
    def wetted_area(self) -> float:
        wetted_area = 0

        if self.nose_cone:
            wetted_area += self.nose_cone.wetted_area
        if self.fins:
            wetted_area += self.fins.wetted_area

        return wetted_area

    @property
    def static_center_of_pressure(self) -> float:
        """
        Calculate the static center of pressure of the rocket.
        Implementation derived from the Barrowman equations.
        """
        if not self.nose_cone or not self.fins:
            raise ValueError("Nose cone and fins must be defined to calculate the static center of pressure")

        # Rocket geometry
        L_N = self.nose_cone.length
        d = self.nose_cone.diameter
        C_R = self.fins.root_chord
        C_T = self.fins.tip_chord
        S = self.fins.semi_span
        X_R = self.fins.sweep_length
        L_F = sqrt(S**2 + (C_T/2 - C_R/2 + X_R)**2)
        R = self.diameter / 2
        X_B = self.fins.position + L_N
        N = self.fins.num_fins

        # Nose cone terms
        C_NN = 2
        X_N = NoseCone.X_N

        # Fin terms
        C_NF = (1 + R/(S + R)) * (4*N*(S/d)**2/(1 + sqrt(1 + (2*L_F/(C_R + C_T))**2)))
        X_F = X_B + X_R/3 * (C_R + 2*C_T)/(C_R + C_T) + 1/6*(C_R + C_T - C_R*C_T/(C_R + C_T))

        # Sum up coefficients
        C_NR = C_NN + C_NF

        # CP distance from the nose cone tip
        return (C_NN*X_N + C_NF*X_F)/C_NR


    def add_motor(self, position, mass, peak_thrust, burn_time) -> None:
        self.motor = Motor(mass, position, peak_thrust, burn_time)

    def add_motor_from_thrust_curve(self, position, filepath) -> None:
        self.motor = Motor.from_thrust_curve(position, filepath)

    def add_nose_cone(self, position, mass, diameter, length) -> None:
        self.nose_cone = NoseCone(position, mass, diameter, length)

    def add_fins(self, position, num_fins, width, root_chord, tip_chord, semi_span, sweep_length) -> None:
        self.fins = Fins(position, num_fins, width, root_chord, tip_chord, semi_span, sweep_length)

    def thrust(self, t) -> float:
        return self.motor.thrust(t)


if __name__ == "__main__":
    pass
