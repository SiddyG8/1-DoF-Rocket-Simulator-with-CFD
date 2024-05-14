import constants as c
from numpy import e, sin, log10
from math import sqrt


def calculate_mass(total_mass, fuel_mass, delta_mass, t) -> float:
    return max(total_mass - fuel_mass, total_mass + delta_mass * t)


def calculate_gravitational_acceleration(h) -> float:
    return c.GM / pow((c.RE + h), 2)


def calculate_gravitational_force(h, m) -> float:
    g = calculate_gravitational_acceleration(h)
    return -g * m


def calculate_mach_number(v, air_temperature) -> float:
    return abs(v) / sqrt(c.gamma * c.R * air_temperature)


def calculate_air_temperature(h) -> float:
    return c.T0 - c.L * h


def calculate_air_density(air_pressure, air_temperature) -> float:
    return air_pressure / (c.R * air_temperature)


def calculate_air_pressure(air_temperature) -> float:
    return c.p0 * pow((air_temperature / c.T0), (-c.g0 / (c.R * c.L)))


def calculate_drag_coefficient(M) -> float:
    return pow(e, -1.2 * M) * sin(M) + (M / 6) * log10(M + 1)


def calculate_drag_force(air_density, v, A, cd) -> float:
    drag = 0.5 * air_density * v * v * A * cd
    if v >= 0:
        return -drag

    return drag


def within_tolerance(value, target, tolerance):
    return target - tolerance <= value <= target + tolerance
