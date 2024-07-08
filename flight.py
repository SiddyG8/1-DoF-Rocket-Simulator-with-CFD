from rocket import Rocket
import constants as c
import functions as f
import matplotlib.pyplot as plt
import csv


"""
Represents the flight of a rocket.

Attributes:
    rocket (Rocket): The rocket object.
    dt (float): The time step for the simulation.
    t (list[float]): The time values.
    z (list[float]): The altitude values.
    vz (list[float]): The vertical velocity values.
    az (list[float]): The vertical acceleration values.
    masses (list[float]): The mass values.
    M (list[float]): The Mach number values.
    gravitational_forces (list[float]): The gravitational force values.
    drag_forces (list[float]): The drag force values.
    drag_coefficients (list[float]): The drag coefficient values.
    air_pressures (list[float]): The air pressure values.
    air_densities (list[float]): The air density values.
    air_temperatures (list[float]): The air temperature values.
    thrust_forces (list[float]): The thrust force values.
    event_log (dict[str, list[float, str]]): The event log.
    twr (list[float]): The thrust-to-weight ratio values.
    dynamic_pressures (list[float]): The dynamic pressure values.
    mach_to_cd (dict[float, float]): The Mach number to drag coefficient mapping.
"""
class Flight:
    def __init__(self, rocket: Rocket, *, drag_data: str = "", dt: float = 0.05) -> None:
        """
        Initialises the flight simulation.

        Parameters:
            rocket (Rocket): The rocket object.
            drag_data (str): The path to the drag data file.
            dt (float): The time step for the simulation.
        """
        # Initialise parameters
        self.rocket = rocket
        self.dt = dt
        self.t = [0]
        self.z = [0]
        self.vz = [0]
        self.az = [self.rocket.thrust(0) / self.rocket.total_mass]
        self.masses = [self.rocket.total_mass]
        self.M = [f.calculate_mach_number(0, c.T0)]
        self.gravitational_forces = [-c.g0 * self.rocket.total_mass]
        self.drag_forces = [0]
        self.drag_coefficients = [0]
        self.air_pressures = [c.p0]
        self.air_densities = [c.rho0]
        self.air_temperatures = [c.T0]
        self.thrust_forces = [self.rocket.thrust(0)]
        self.event_log = {
            "Motor Ignition": [0, "green"],
            "Motor Burnout": [self.rocket.motor.burn_time, "orange"],
            "Apogee": [0, "blue"],
            "Ground Hit": [0, "red"]
        }
        self.twr = [self.rocket.thrust(0) / (self.rocket.total_mass * c.g0)]
        self.dynamic_pressures = [0]
        self.mach_to_cd = {}

        # Initialise the CFD drag data if provided
        if drag_data:
            self.initialise_drag_data(drag_data)

        # Simulate the flight based off initial parameters
        self.simulate()

    def initialise_drag_data(self, drag_data: str) -> None:
        """
        Initialises the drag data from a CSV file.

        Parameters:
            drag_data (str): The path to the drag data file.
        """
        with open(drag_data) as file:
            reader = csv.DictReader(file)
            for row in reader:
                pitch, roll, flap = row["Pitch (deg)"], row["Roll (deg)"], row["Flap (% extension)"]
                if pitch == 0 and roll == 0 and flap == 0:
                    M, cd = row["Mach Number"], row["Total C_D"]
                    self.mach_to_cd[float(M)] = float(cd)

    def f(self, t: float, s: float, v: float) -> float:
        """
        Calculates the net acceleration of the rocket.

        Parameters:
            t (float): The time.
            s (float): The altitude.
            v (float): The velocity.

        Returns:
            float: The net acceleration of the rocket.
        """
        # Calculate parameters required for the force calculations
        air_temperature = f.calculate_air_temperature(s)
        air_pressure = f.calculate_air_pressure(air_temperature)
        air_density = f.calculate_air_density(air_pressure, air_temperature)
        mach_number = f.calculate_mach_number(v, air_temperature)
        drag_coefficient = f.calculate_drag_coefficient(mach_number, self.mach_to_cd)

        # Calculate the net force acting on the rocket
        net_force = \
            self.rocket.thrust(t) \
            + f.calculate_drag_force(air_density, v, self.rocket.wetted_area, drag_coefficient) \
            + f.calculate_gravitational_force(s, self.masses[-1])

        # Return the net acceleration
        return net_force / self.masses[-1]

    def simulate(self) -> None:
        """
        Simulates the flight of the rocket.
        """
        # Simulation loop
        apogee_reached = False
        while self.z[-1] >= 0:
            # Calculate the next state
            t = self.t[-1] + self.dt
            z = self.z[-1] + self.dt * self.vz[-1]
            vz = self.vz[-1] + self.dt * self.az[-1]
            az = self.f(t, z, vz)
            mass = f.calculate_mass(self.rocket.total_mass, self.rocket.motor.mass, self.rocket.motor.delta_mass, t)
            gravitational_force = f.calculate_gravitational_force(z, mass)
            air_temperature = f.calculate_air_temperature(z)
            air_pressure = f.calculate_air_pressure(air_temperature)
            air_density = f.calculate_air_density(air_pressure, air_temperature)
            M = f.calculate_mach_number(vz, air_temperature)
            cd = f.calculate_drag_coefficient(M)
            drag_force = f.calculate_drag_force(air_density, vz, self.rocket.wetted_area, cd)
            thrust_force = self.rocket.thrust(t)
            twr = thrust_force / (mass * c.g0)
            dynamic_pressure = f.calculate_dynamic_pressure(air_density, vz)

            # Apogee event
            if vz <= 0 and not apogee_reached:
                self.event_log["Apogee"][0] = t
                apogee_reached = True

            # Update parameters
            self.t.append(t)
            self.z.append(z)
            self.vz.append(vz)
            self.az.append(az)
            self.masses.append(mass)
            self.gravitational_forces.append(gravitational_force)
            self.air_temperatures.append(air_temperature)
            self.air_densities.append(air_density)
            self.air_pressures.append(air_pressure)
            self.M.append(M)
            self.drag_coefficients.append(cd)
            self.drag_forces.append(drag_force)
            self.thrust_forces.append(thrust_force)
            self.twr.append(twr)
            self.dynamic_pressures.append(dynamic_pressure)

        # Ground hit event
        self.event_log["Ground Hit"][0] = self.t[-1]
    
    def off_rod_velocity(self, rod_length: float) -> tuple[float, float]:
        """
        Calculates the velocity of the rocket when it leaves the launch rod.

        Parameters:
            rod_length (float): The length of the launch rod.

        Returns:
            tuple[float, float]: The time and velocity when the rocket leaves the launch rod.
        """
        for i in range(len(self.z)):
            if self.z[i] >= rod_length:
                return self.t[i], self.vz[i]

    def plot(self, x: str, *y: tuple[str, str] | str, x_label: str = "x", y_label: str = "y", events=True) -> None:
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

        # Label the axes
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # Plot events if required
        if events:
            # Plot events as vertical lines with labels and timestamps
            for event, details in self.event_log.items():
                t, color = details
                plt.axvline(t, color=color, linestyle="--", linewidth=1)
                plt.text(t, int(sum(plt.gca().get_ylim()[:2]) / 2), event, color=color, fontsize=8, rotation=-90, verticalalignment="center")
                plt.text(t, plt.gca().get_ylim()[1], f"{t:.2f}s", color=color, fontsize=8, ha="center", va="bottom")

        # Show legend if required
        if legend:
            plt.legend()

        # Show the plot with gridlines
        plt.grid()
        plt.show()


if __name__ == "__main__":
    pass
