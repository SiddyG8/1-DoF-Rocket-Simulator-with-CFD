import constants as c
from math import sqrt, e, sin, log10


def calculate_mass(total_mass, fuel_mass, delta_mass, t) -> float:
    """
    Calculates the mass of the rocket at time t.

    Parameters:
        total_mass (float): The total mass of the rocket.
        fuel_mass (float): The mass of the fuel.
        delta_mass (float): The mass flow rate of the fuel.
        t (float): The time.

    Returns:
        float: The mass of the rocket at time t.
    """
    return max(total_mass - fuel_mass, total_mass - delta_mass * t)


def calculate_gravitational_acceleration(h) -> float:
    """
    Calculates the gravitational acceleration at a given altitude.

    Parameters:
        h (float): The altitude.

    Returns:
        float: The gravitational acceleration at the given altitude.
    """
    return c.GM / pow((c.RE + h), 2)


def calculate_gravitational_force(h, m) -> float:
    """
    Calculates the gravitational force acting on an object at a given altitude.

    Parameters:
        h (float): The altitude.
        m (float): The mass of the object.

    Returns:
        float: The gravitational force acting on the object at the given altitude
    """
    g = calculate_gravitational_acceleration(h)
    return -g * m


def calculate_mach_number(v, air_temperature) -> float:
    """
    Calculates the Mach number of an object moving at a given velocity.

    Parameters:
        v (float): The velocity of the object.
        air_temperature (float): The temperature of the air.

    Returns:
        float: The Mach number of the object.
    """
    return abs(v) / sqrt(c.gamma * c.R * air_temperature)


def calculate_air_temperature(h) -> float:
    """
    Calculates the temperature of the air at a given altitude.

    Parameters:
        h (float): The altitude.

    Returns:
        float: The temperature of the air at the given altitude.
    """
    return c.T0 - c.L * h


def calculate_air_density(air_pressure, air_temperature) -> float:
    """
    Calculates the density of the air at a given pressure and temperature.

    Parameters:
        air_pressure (float): The pressure of the air.
        air_temperature (float): The temperature of the air.

    Returns:
        float: The density of the air.
    """
    return air_pressure / (c.R * air_temperature)


def calculate_air_pressure(air_temperature) -> float:
    """
    Calculates the pressure of the air at a given temperature.

    Parameters:
        air_temperature (float): The temperature of the air.

    Returns:
        float: The pressure of the air
    """
    return c.p0 * pow((air_temperature / c.T0), (-c.g0 / (c.R * c.L)))


def calculate_dynamic_pressure(air_density, u) -> float:
    """
    Calculates the dynamic pressure of the air.

    Parameters:
        air_density (float): The density of the air.
        u (float): The velocity of the object.

    Returns:
        float: The dynamic pressure
    """
    return 0.5 * air_density * u * u


def calculate_drag_coefficient(M: float, drag_data: dict[float, float] = None) -> float:
    """
    Calculates the drag coefficient of an object moving at a given Mach number.

    Parameters:
        M (float): The Mach number of the object.
        drag_data (dict[float, float]): The drag data.

    Returns:
        float: The drag coefficient
    """
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
    """
    Calculates the drag force acting on an object moving at a given velocity.

    Parameters:
        air_density (float): The density of the air.
        v (float): The velocity of the object.
        A (float): The reference area of the object.
        cd (float): The drag coefficient of the object.

    Returns:
        float: The drag force acting on the object.
    """
    drag = 0.5 * air_density * v * v * A * cd

    # Return negative drag if the velocity is negative
    return drag if v < 0 else -drag
