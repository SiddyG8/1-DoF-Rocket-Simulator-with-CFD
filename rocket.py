from dataclasses import dataclass
from math import pi, e, log, sqrt
from constants import g0


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


class Parachute:
    def __init__(self, diameter: float, drag_coefficient: float, mass: float) -> None:
        self.diameter = diameter
        self.drag_coefficient = drag_coefficient
        self.mass = mass

    @property
    def wetted_area(self) -> float:
        return pi * (self.diameter / 2) ** 2


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
        self.main_parachute = None
        self.drogue_parachute = None

    @property
    def total_mass(self) -> float:
        total_mass = self.body_mass + self.motor.mass

        if self.nose_cone:
            total_mass += self.nose_cone.mass
        if self.fins:
            total_mass += self.fins.total_mass
        if self.main_parachute:
            total_mass += self.main_parachute.mass
        if self.drogue_parachute:
            total_mass += self.drogue_parachute.mass

        return total_mass

    @property
    def dry_mass(self) -> float:
        return self.total_mass - self.motor.mass

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
        X_B = self.body_length + L_N - self.fins.position
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

    @property
    def initial_center_of_gravity(self) -> float:
        if not self.nose_cone or not self.fins:
            raise ValueError("Nose cone and fins must be defined to calculate the static center of pressure")

        X_N = self.nose_cone.length * self.nose_cone.X_N + self.body_length
        X_F = self.fins.position
        X_R = self.body_length / 2
        X_M = self.motor.length / 2

        # fuel_mass = rocket.mass - current_mass

        X_CG = (X_N * self.nose_cone.mass + X_F * self.fins.total_mass + X_R * self.total_mass + X_M * self.motor.mass) / self.total_mass

        return self.body_length + self.nose_cone.length - X_CG

    @property
    def static_margin(self) -> float:
        return (self.static_center_of_pressure - self.initial_center_of_gravity) / self.diameter

    @property
    def initial_twr(self) -> float:
        return self.motor.thrust(0) / (self.total_mass * g0)

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
    
    def attach_main_parachute(self, parachute: Parachute) -> None:
        """
        Attaches a main parachute to the rocket.
        """
        self.main_parachute = parachute

    def attach_drogue_parachute(self, parachute: Parachute) -> None:
        """
        Attaches a drogue parachute to the rocket.
        """
        self.drogue_parachute = parachute

    def thrust(self, t) -> float:
        """
        Returns the thrust of the rocket at a given time. Delegates to the motor's thrust method.
        """
        return self.motor.thrust(t)

    def __str__(self) -> str:
        outstring = ""
        outstring += f"Mass: {self.total_mass:.3f} kg"
        outstring += f"\nDiameter: {self.diameter:.3f} m"
        outstring += f"\nLength: {self.body_length:.3f} m"
        outstring += f"\nInitial CP (from nose): {self.static_center_of_pressure:.3f} m"
        outstring += f"\nInitial CG (from nose): {self.initial_center_of_gravity:.3f} m"
        outstring += f"\nStatic Margin: {self.static_margin:.2f}"
        outstring += f"\nInitial TWR: {self.initial_twr:.2f}"
        return outstring


if __name__ == "__main__":
    pass
