from rocket import Rocket
import constants as c
from numpy import e, sin, log10

class Flight:
    """
    A class to simulate the flight of a rocket.

    Attributes
    ----------
    rocket : Rocket
        The rocket to simulate.
    dt : float
        The time step of the simulation (s).
    t : float
        The current time of the simulation (s).
    z : float
        The current altitude of the rocket (m).
    vz : float
        The current vertical velocity of the rocket (m/s).
    az : float
        The current vertical acceleration of the rocket (m/s^2).
    phase : str
        The current phase of the flight.
    """
    def __init__(
        self,
        rocket: Rocket,
        dt: float
    ) -> None:
        """
        Initialize the flight simulation.

        Parameters
        ----------
        rocket : Rocket
            The rocket to simulate.
        dt : float
            The time step of the simulation (s).

        Returns
        -------
        None
        """
        self.rocket = rocket
        self.dt = dt
        self.t = 0
        self.z = 0
        self.vz = 0
        self.az = 0
        self.phase = "pad"

    @property
    def altitude(self) -> float:
        return self.z

    @property
    def drag_coefficient(self) -> float:
        return pow(e, -1.2 * self.M) * sin(self.M) + (self.M / 6) * log10(self.M + 1)

    @property
    def mach_number(self) -> float:
        pass

    def simulate(self) -> None:
        """
        Simulate the flight of the rocket.
        """
        pass

if __name__ == "__main__":
    pass
