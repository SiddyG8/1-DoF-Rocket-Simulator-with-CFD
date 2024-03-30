from numpy import e, sin, log10

class Functions:
    @classmethod
    def drag_coefficient(M: float) -> float:
        """
        Calculate the drag coefficient of the rocket.

        Parameters:
            M (float): Mach number of the rocket.
        Returns:
            float: Drag coefficient of the rocket.
        """
        return pow(e, -1.2 * M) * sin(M) + (M / 6) * log10(M + 1)

    @classmethod
    def drag_force(Cd: float, S: float, p: float, v: float) -> float:
        """
        Calculate the drag force of the rocket.

        Parameters:
            Cd (float): Drag coefficient of the rocket.
            S (float): Wetted area of the rocket.
            p (float): Air density.
            v (float): Velocity of the rocket.
        Returns:
            float: Drag force on the rocket.
        """
        return 0.5 * p * v * v * S * Cd
    
    @classmethod
    def wetted_area(D: float, L: float) -> float:
        """
        Calculate the wetted area of the rocket.

        Parameters:
            D (float): Diameter of the rocket.
            L (float): Length of the rocket.
        Returns:
            float: Wetted area of the rocket.
        """
        return D * L

if __name__ == "__main__":
    pass
