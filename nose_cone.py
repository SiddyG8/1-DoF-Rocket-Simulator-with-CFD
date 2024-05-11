import math

class NoseCone:
    def __init__(self, length, diameter) -> None:
        self.length = length
        self.diameter = diameter

    @property
    def wetted_area(self) -> float:
        return math.pi * (self.diameter / 2) ** 2

if __name__ == "__main__":
    pass
