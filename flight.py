from rocket import Rocket
import constants as c
from math import sqrt, e, sin, log10
import matplotlib.pyplot as plt
import numpy as np
import csv
import os


class Flight:
    def __init__(self, rocket: Rocket, *, rod_height: float = 3, drag_data: str = "", dt: float = 0.05) -> None:
        # Flight parameters
        self.rocket = rocket
        self.rod_height = rod_height
        self.drag_data = drag_data
        self.dt = dt
        self.iterations = 1

        # Initial conditions
        self.t = 0
        self.z = 0
        self.vz = 0
        self.az = self.rocket.thrust(self.t) / self.rocket.total_mass - c.g0
        self.mass = self.rocket.total_mass
        self.air_temperature = c.T0
        self.air_pressure = c.p0
        self.air_density = c.rho0
        self.M = 0
        self.g = -c.g0
        self.gravitational_force = -c.g0 * self.rocket.total_mass
        self.cd = 0
        self.drag_force = 0
        self.thrust_force = self.rocket.thrust(self.t)
        self.dynamic_pressure = 0
        self.twr = self.rocket.thrust(self.t) / (self.rocket.total_mass * c.g0)

        # Initialise data lists
        self.times = [self.t]
        self.altitudes = [self.z]
        self.velocities = [self.vz]
        self.accelerations = [self.az]
        self.masses = [self.mass]
        self.gravitational_accelerations = [self.g]
        self.gravitational_forces = [self.gravitational_force]
        self.air_temperatures = [self.air_temperature]
        self.air_pressures = [self.air_pressure]
        self.air_densities = [self.air_density]
        self.mach_numbers = [self.M]
        self.drag_coefficients = [self.cd]
        self.drag_forces = [self.drag_force]
        self.thrust_forces = [self.thrust_force]
        self.dynamic_pressures = [self.dynamic_pressure]
        self.twrs = [self.twr]

        # Event log (populated during simulation)
        self.event_log = {}
        self.flight_stats = {}

        # Initialise drag data if provided
        self.using_drag_data = os.path.isfile(self.drag_data)
        if self.using_drag_data:
            self.load_drag_data()

        # Simulate the flight based off initial parameters
        self.simulate()

    def simulate(self) -> None:
        while self.z >= 0:
            self.update_parameters()
            self.update_events()
            self.log_parameters()

    def load_drag_data(self) -> None:
        with open(self.drag_data) as file:
            reader = csv.DictReader(file)
            self.mach_data = []
            self.cd_data = []
            for row in reader:
                self.mach_data.append(float(row["Mach Number"]))
                self.cd_data.append(float(row["Total C_D"]))

    def update_events(self) -> None:
        # Update event log
        if self.z < self.altitudes[-1] and "Apogee" not in self.event_log:
            self.event_log["Apogee"] = self.t, self.altitudes[-1]
        if self.thrust_force > 0 and "Motor Ignition" not in self.event_log:
            self.event_log["Motor Ignition"] = self.times[-1]
        if self.thrust_force == 0 and "Motor Burnout" not in self.event_log:
            self.event_log["Motor Burnout"] = self.t
        if self.z <= 0 and self.t > 0 and "Ground Hit" not in self.event_log:
            self.event_log["Ground Hit"] = self.t, self.vz
        # Chute deployment events soon

        # Update flight stats
        if self.vz < self.velocities[-1] and "Max Velocity" not in self.flight_stats:
            self.flight_stats["Max Velocity"] = self.t, self.velocities[-1]
        if self.az < self.accelerations[-1] and "Max Acceleration" not in self.flight_stats:
            self.flight_stats["Max Acceleration"] = self.t, self.accelerations[-1]
        if self.dynamic_pressure < self.dynamic_pressures[-1] and "Max Q" not in self.flight_stats:
            self.flight_stats["Max Q"] = self.t, self.dynamic_pressures[-1]
        if self.z >= self.rod_height and "Off-Rod Velocity" not in self.flight_stats:
            self.flight_stats["Off-Rod Velocity"] = self.t, self.vz

    def log_parameters(self) -> None:
        self.times.append(self.t)
        self.altitudes.append(self.z)
        self.velocities.append(self.vz)
        self.accelerations.append(self.az)
        self.masses.append(self.mass)
        self.gravitational_accelerations.append(self.g)
        self.gravitational_forces.append(self.gravitational_force)
        self.air_temperatures.append(self.air_temperature)
        self.air_pressures.append(self.air_pressure)
        self.air_densities.append(self.air_density)
        self.mach_numbers.append(self.M)
        self.drag_coefficients.append(self.cd)
        self.drag_forces.append(self.drag_force)
        self.thrust_forces.append(self.thrust_force)
        self.dynamic_pressures.append(self.dynamic_pressure)
        self.twrs.append(self.twr)

    def update_parameters(self) -> None:
        self.euler_integration()
        self.update_mass()
        self.update_gravitational_acceleration()
        self.update_gravitational_force()
        self.update_air_temperature()
        self.update_air_pressure()
        self.update_air_density()
        self.update_mach_number()
        self.update_drag_coefficient()
        self.update_drag_force()
        self.update_thrust_force()
        self.update_dynamic_pressure()
        self.update_twr()
        self.iterations += 1

    def f(self) -> float:
        return (self.thrust_force + self.drag_force + self.gravitational_force) / self.mass

    def euler_integration(self) -> None:
        self.t += self.dt
        self.az = self.f()
        self.vz += self.az * self.dt
        self.z += self.vz * self.dt

    def update_mass(self) -> None:
        self.mass = max(self.rocket.dry_mass, self.rocket.total_mass - self.rocket.motor.delta_mass * self.t)

    def update_gravitational_acceleration(self) -> None:
        self.g = -c.GM / (c.RE + self.z)**2

    def update_gravitational_force(self) -> None:
        self.gravitational_force = self.mass * self.g

    def update_air_temperature(self) -> None:
        self.air_temperature = c.T0 + c.L * self.z

    def update_air_pressure(self) -> None:
        self.air_pressure = c.p0 * pow(self.air_temperature / c.T0, -c.g0 / (c.R * c.L))

    def update_air_density(self) -> None:
        self.air_density = self.air_pressure / (c.R * self.air_temperature)

    def update_mach_number(self) -> None:
        self.M = abs(self.vz) / sqrt(c.gamma * c.R * self.air_temperature)

    def update_drag_coefficient(self) -> None:
        if self.using_drag_data:
            self.cd = np.interp(self.M, self.mach_data, self.cd_data)
            return None

        self.cd = pow(e, -1.2 * self.M) * sin(self.M) + (self.M / 6) * log10(self.M + 1)

    def update_drag_force(self) -> None:
        drag_force = 0.5 * self.air_density * self.vz**2 * self.rocket.wetted_area * self.cd
        self.drag_force = -drag_force if self.vz >= 0 else drag_force

    def update_thrust_force(self) -> None:
        self.thrust_force = self.rocket.thrust(self.t)

    def update_dynamic_pressure(self) -> None:
        self.dynamic_pressure = 0.5 * self.air_density * self.vz**2

    def update_twr(self) -> None:
        self.twr = abs(self.thrust_force / (self.mass * c.g0))

    def plot(self, x: str, *y: tuple[str, str] | str, title: str = "Flight Data", x_label: str = "x", y_label: str = "y", events=True) -> None:
        """
        Plots the simulation data.

        Parameters:
            x (str): The x-axis data.
            y (tuple[str, str] | str): The y-axis data.
            x_label (str): The x-axis label.
            y_label (str): The y-axis label.
            events (bool): Whether to plot events.

        Example:
            flight.plot(
                "t",
                ("z", "Altitude"),
                ("vz", "Velocity (m/s)"),
                ("az", "Acceleration (m/s^2)"),
                x_label="Time (s)",
                y_label="Altitude (m)"
            )
            This will plot the altitude, velocity, and acceleration against time. It will also label the curves and axes.
        """
        # Get the x-axis data via its attribute name
        x_var = getattr(self, x)

        # Plot the data and label the curves
        legend = False
        for var in y:
            y_attr, label = var if isinstance(var, tuple) else (var, None)
            legend = True if label else legend
            plt.plot(x_var, getattr(self, y_attr), label=label)

        # Label the axes and title
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # Plot events if required
        if events:
            # Plot events as vertical lines with labels and timestamps
            for event, value in self.event_log.items():
                t, _ = value if isinstance(value, tuple) else (value, None)
                plt.axvline(t, linestyle="--", linewidth=1)
                plt.text(t, int(sum(plt.gca().get_ylim()[:2]) / 2), event, fontsize=8, rotation=-90, verticalalignment="center")
                plt.text(t, plt.gca().get_ylim()[1], f"{t:.2f}s", fontsize=8, ha="center", va="bottom")

        # Show legend if required
        if legend:
            plt.legend()

        # Show the plot with gridlines
        plt.grid()
        plt.show()

    def write_to_file(self, output_path: str) -> None:
        data_dict = {
            "Time (s)": self.times,
            "Altitude (m)": self.altitudes,
            "Velocity (m/s)": self.velocities,
            "Acceleration (m/s^2)": self.accelerations,
            "Mass (kg)": self.masses,
            "Gravitation (m/s^2)": self.gravitational_accelerations,
            "Gravitational Force (N)": self.gravitational_forces,
            "Air Temperature (K)": self.air_temperatures,
            "Air Pressure (Pa)": self.air_pressures,
            "Air Density (kg/m^3)": self.air_densities,
            "Mach Number": self.mach_numbers,
            "Drag Coefficient": self.drag_coefficients,
            "Drag Force (N)": self.drag_forces,
            "Thrust Force (N)": self.thrust_forces,
            "Dynamic Pressure (Pa)": self.dynamic_pressures,
            "TWR": self.twrs
        }

        with open(output_path, "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=data_dict.keys())

            # Write the header to the file
            writer.writeheader()

            # Write the parameters to the file
            for i in range(self.iterations):
                writer.writerow({key: data[i] for key, data in data_dict.items()})

    def summary(self) -> None:
        print("Flight Summary")
        print("---------------------------------")
        print(f"Apogee: {self.event_log["Apogee"][1]:.2f}m at {self.event_log["Apogee"][0]:.2f}s")
        print(f"Motor Ignition: {self.event_log["Motor Ignition"]:.2f}s")
        print(f"Motor Burnout: {self.event_log["Motor Burnout"]:.2f}s")
        print(f"Flight Time: {self.event_log["Ground Hit"][0]:.2f}s with a landing velocity of {self.event_log["Ground Hit"][1]:.2f}m/s")
        print(f"Max Velocity: {self.flight_stats["Max Velocity"][1]:.2f}m/s at {self.flight_stats["Max Velocity"][0]:.2f}s")
        print(f"Max Acceleration: {self.flight_stats["Max Acceleration"][1]:.2f}m/s^2 at {self.flight_stats["Max Acceleration"][0]:.2f}s")
        print(f"Max Q: {self.flight_stats["Max Q"][1]/1000:.2f}kPa at {self.flight_stats["Max Q"][0]:.2f}s")
        print(f"Off-Rod Velocity: {self.flight_stats["Off-Rod Velocity"][1]:.2f}m/s at {self.flight_stats["Off-Rod Velocity"][0]:.2f}s")
        print("---------------------------------")
        print()


if __name__ == "__main__":
    pass
