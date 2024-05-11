from rocket import Rocket
import constants as c
from integerator import euler
from numpy import e, sin, log10
import math

class Flight:
    def __init__(self, rocket: Rocket, dt: float = 0.05) -> None:
        self.rocket = rocket
        self.dt = dt
        self.step = 0
        self.t = [0]
        self.z = [0]
        self.vz = [0]
        self.az = [0]
        self.M = [0]
        self.g = [c.g0]
        self.drag_forces = [0]
        self.drag_coefficients = [0]
        self.air_pressures = [c.p0]
        self.air_densities = [c.rho0]
        self.air_temperatures = [c.T0]
        self.thrust_forces = [0]
        self.phases = ["on pad"]

        # Run the simulation
        self.simulate()

    def calculate_gravitational_acceleration(self, step) -> float:
        return c.GM / (c.RE + self.z[step]) ** 2

    def calculate_gravitational_force(self, step) -> float:
        return -self.rocket.mass * self.g[step]

    def calculate_thrust_force(self, step) -> float:
        return self.rocket.thrust(self.t[step])

    def calculate_mach_number(self, step) -> float:
        return abs(self.vz[step]) / math.sqrt(c.gamma * c.R * self.air_temperatures[step])

    def calculate_altitude(self, step) -> float:
        return self.z[step]

    def calculate_air_temperature(self, step) -> float:
        return c.T0 - c.L * self.z[step]

    def calculate_air_density(self, step) -> float:
        return self.air_pressures[step] / (c.R * self.air_temperatures[step])

    def calculate_pressure(self, step) -> float:
        return c.p0 * (self.air_temperatures[step] / c.T0) ** (-c.g0 / (c.R * c.L))

    def calculate_drag_coefficient(self, step) -> float:
        return pow(e, -1.2 * self.M[step]) * sin(self.M[step]) + (self.M[step] / 6) * log10(self.M[step] + 1)

    def calculate_drag_force(self, step) -> float:
        return 0.5 * self.air_densities[step] * self.vz[step] * self.vz[step] * self.rocket.wetted_area * self.drag_coefficients[step]
    
    def f(self, t, s, v):
        self.step += 1
        self.t.append(t)
        self.z.append(s)
        self.vz.append(v)
        return self.calculate_thrust_force(self.step) - self.calculate_drag_force(self.step) - self.calculate_gravitational_force(self.step)

    def simulate(self) -> None:
        while self.phases[self.step] != "landed":
            euler(self.f, self.dt)

    def plot(*variables, x_label, y_label) -> None:
        pass

if __name__ == "__main__":
    pass
