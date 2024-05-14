from math import e, log

class Motor:
    def __init__(self, T0: float, tb: float) -> None:
        self.T0 = T0
        self.tb = tb

    @classmethod
    def from_thrust_curve(cls, thrust_curve_path: str) -> None:
        raise NotImplementedError

    def thrust(self, t) -> float:
        return max(0, self.T0 * (1 - (10**-5) * (e**(log(10**5)/self.tb * t))))
