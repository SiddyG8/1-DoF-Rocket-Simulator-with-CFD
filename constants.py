import numpy
from dataclasses import dataclass

@dataclass
class Constants:
    g0: float = 9.80665, # Mean gravitational acceleration at sea level (m/s^2)
    R: float = 287.16, # Specfic gas constant for air
    p0: int = 101325, # Standard sea level pressure (Pa)
    T0: float = 288.15, # Standard sea level temperature (K)
    L: float = -0.0065, # Lapse rate (K/m)

