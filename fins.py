class Fins:
    def __init__(self, root_chord, tip_chord, span, sweep_angle, thickness, no_fins) -> None:
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.span = span
        self.sweep_angle = sweep_angle
        self.thickness = thickness
        self.no_fins = no_fins

    @property
    def wetted_area(self) -> float:
        return self.span * self.thickness
