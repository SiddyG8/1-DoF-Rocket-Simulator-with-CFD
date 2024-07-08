import constants as c
from math import sqrt, e, sin, log10


def calculate_mass(total_mass, fuel_mass, delta_mass, t) -> float:
    return max(total_mass - fuel_mass, total_mass - delta_mass * t)


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


def calculate_dynamic_pressure(air_density, u) -> float:
    return 0.5 * air_density * u * u


def calculate_drag_coefficient(M: float, drag_data: dict[float, float] = None) -> float:
    if drag_data:
        for data_M, data_cd in drag_data.items():
            if M - data_M > 0:
                M1 = data_M
                cd1 = data_cd
            elif M - data_M < 0:
                M2 = data_M
                cd2 = data_cd
            else:
                return data_cd

        return (cd1*(M2 - M) + cd2*(M - M1))/(M2 - M1)

    return pow(e, -1.2 * M) * sin(M) + (M / 6) * log10(M + 1)


def calculate_drag_force(air_density, v, A, cd) -> float:
    drag = 0.5 * air_density * v * v * A * cd

    # Return negative drag if the velocity is negative
    return drag if v < 0 else -drag
