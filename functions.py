from numpy import e, sin, log10

class Functions:
    def drag_coefficient(
        M: float, # Mach number
    ) -> float:
        return pow(e, -1.2 * M) * sin(M) + (M / 6) * log10(M + 1)

if __name__ == "__main__":
    pass
