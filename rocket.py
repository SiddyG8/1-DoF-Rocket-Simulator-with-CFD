from dataclasses import dataclass
from math import pi, e, log


"""
Represents a rocket motor.

Attributes:
    mass (float): The mass of the motor.
    length (float): The length of the motor.
    peak_thrust (float): The peak thrust of the motor.
    burn_time (float): The burn time of the motor.
"""
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
        """
        Returns the thrust of the motor at a given time. Equation is based on the motor's peak thrust and burn time.
        """
        # Equation driven thrust curve that estimates a motor's thrust over time based on its peak thrust and burn time
        return max(0, self.peak_thrust * (1 - (10**-5) * (e**(log(10**5)/self.burn_time * t))))


"""
Represents the (Von Kármán) nose cone of a rocket.

Attributes:
    mass (float): The mass of the nose
    diameter (float): The diameter of the nose
    length (float): The length of the nose
"""
class NoseCone:
    X_N = 0.466  # Nose cone length ratio (Ogive)

    def __init__(self, mass: float, diameter: float, length: float) -> None:
        """
        Initialises the nose cone.

        Parameters:
            mass (float): The mass of the nose cone.
            diameter (float): The diameter of the nose cone.
            length (float): The length of the nose cone.
        """
        self.mass = mass
        self.diameter = diameter
        self.length = length

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2


"""
Represents the fins of a rocket.

Attributes:
    position (float): The position of the fins.
    fin_mass (float): The mass of the fins.
    num_fins (int): The number of fins.
    width (float): The width of the fins.
    root_chord (float): The root chord of the fins.
    tip_chord (float): The tip chord of the fins.
    semi_span (float): The semi span of the fins.
    sweep_length (float): The sweep length of the fins.
"""
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


"""
Represents a rocket.

Attributes:
    body_mass (float): The mass of the rocket body.
    diameter (float): The diameter of the rocket body.
    body_length (float): The length of the rocket body.
    motor (Motor): The motor of the rocket.
    nose_cone (NoseCone): The nose cone of the rocket.
    fins (Fins): The fins of the rocket.
"""
class Rocket:
    def __init__(self, motor: Motor, body_mass: float, diameter: float, body_length: float) -> None:
        """
        Initialises the rocket.

        Parameters:
            motor (Motor): The motor of the rocket.
            body_mass (float): The mass of the rocket body.
            diameter (float): The diameter of the rocket body.
            body_length (float): The length of the rocket
        """
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
        """
        Attaches a nose cone to the rocket.
        """
        self.nose_cone = nose_cone

    def attach_fins(self, fins: Fins) -> None:
        """
        Attaches fins to the rocket.
        """
        self.fins = fins

    def thrust(self, t) -> float:
        """
        Returns the thrust of the rocket at a given time. Delegates to the motor's thrust method.
        """
        return self.motor.thrust(t)


if __name__ == "__main__":
    pass
