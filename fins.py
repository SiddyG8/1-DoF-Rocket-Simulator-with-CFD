class Fins:
    def __init__(self, num_fins, width, root_chord, tip_chord, span, sweep_angle) -> None:
        self.num_fins = num_fins
        self.width = width
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.span = span
        self.sweep_angle = sweep_angle

    @property
    def wetted_area(self) -> float:
        return self.num_fins * self.width * self.span
