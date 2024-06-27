from dataclasses import dataclass


@dataclass
class Fins:
    num_fins: int
    width: float
    root_chord: float
    tip_chord: float
    semi_span: float
    sweep_length: float

    @property
    def wetted_area(self) -> float:
        return self.num_fins * self.width * self.semi_span
