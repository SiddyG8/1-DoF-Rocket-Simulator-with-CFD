from typing import ClassVar
from dataclasses import dataclass

@dataclass
class NoseCone:
    # Constants
    C_NN: ClassVar[int] = 2

    # Parameters
    length: float
    diameter: float
    type: str = "ogive"
