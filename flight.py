from rocket import Rocket
import constants as c
import functions as f
import matplotlib.pyplot as plt

class Flight:
    def __init__(self, rocket: Rocket, dt: float = 0.05) -> None:
        self.rocket = rocket
        self.dt = dt
        self.t = [0]
        self.z = [0]
        self.vz = [0]
        self.az = [self.rocket.thrust(0) / self.rocket.mass]
        self.masses = [self.rocket.mass]
        self.M = [f.calculate_mach_number(0, c.T0)]
        self.gravitional_forces = [-c.g0 * self.rocket.mass]
        self.drag_forces = [0]
        self.drag_coefficients = [0]
        self.air_pressures = [c.p0]
        self.air_densities = [c.rho0]
        self.air_temperatures = [c.T0]
        self.thrust_forces = [self.rocket.thrust(0)]
        self.phases = ["on pad"]

        # Simulate the flight based off initial parameters
        self.simulate()

    def f(self, t, s, v):
        air_temperature = f.calculate_air_temperature(s)
        air_pressure = f.calculate_air_pressure(air_temperature)
        air_density = f.calculate_air_density(air_pressure, air_temperature)
        mach_number = f.calculate_mach_number(v, air_temperature)
        drag_coefficient = f.calculate_drag_coefficient(mach_number)
        net_force = self.rocket.thrust(t) + f.calculate_drag_force(air_density, v, self.rocket.wetted_area, drag_coefficient) + f.calculate_gravitational_force(s, self.masses[-1])
        return net_force / self.masses[-1]

    def simulate(self) -> None:
        while "landed" not in self.phases:
            t = self.t[-1] + self.dt
            z = self.z[-1] + self.dt * self.vz[-1]
            vz = self.vz[-1] + self.dt * self.az[-1]
            az = self.f(t, z, vz)
            mass = f.calculate_mass(self.rocket.mass, self.rocket.fuel_mass, self.rocket.delta_mass, t)
            gravitional_force = f.calculate_gravitational_force(z, mass)
            air_temperature = f.calculate_air_temperature(z)
            air_pressure = f.calculate_air_pressure(air_temperature)
            air_density = f.calculate_air_density(air_pressure, air_temperature)
            M = f.calculate_mach_number(vz, air_temperature)
            cd = f.calculate_drag_coefficient(M)
            drag_force = f.calculate_drag_force(air_density, vz, self.rocket.wetted_area, cd)
            thrust_force = self.rocket.thrust(t)

            if z < 0:
                phase = "landed"
            else:
                phase = "on pad"

            self.t.append(t)
            self.z.append(z)
            self.vz.append(vz)
            self.az.append(az)
            self.masses.append(mass)
            self.gravitional_forces.append(gravitional_force)
            self.air_temperatures.append(air_temperature)
            self.air_densities.append(air_density)
            self.air_pressures.append(air_pressure)
            self.M.append(M)
            self.drag_coefficients.append(cd)
            self.drag_forces.append(drag_force)
            self.thrust_forces.append(thrust_force)
            self.phases.append(phase)

    def plot(self, x, *y: tuple[str, str], x_label: str = "x", y_label: str = "y") -> None:
        x_var = getattr(self, x)
        legend_needed = False
        for var in y:
            y_attr, label = var if isinstance(var, tuple) else (var, None)
            legend_needed = True if label else legend_needed
            plt.plot(x_var, getattr(self, y_attr), label=label)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        if legend_needed:
            plt.legend()
        plt.show()
        # TODO: Create a seperate class for handling plots
        # TODO: Create a seperate class for displaying flight statistics


if __name__ == "__main__":
    pass
