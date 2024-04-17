from numpy import e, log

class Motor:
    def __init__(self, T0, tb) -> None:
        self.T0 = T0
        self.tb = tb

    def thrust(self, t) -> float:
        return self.T0 * (1 - (10**-5) * (e**(log(10**5)/self.tb * t)))

if __name__ == "__main__":
    pass
