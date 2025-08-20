from .ligand import Ligand
from abc import ABC, abstractclassmethod

class Substructure:
    def __init__(self) -> None:
        pass

class Pattern(ABC):
    def __init__(self) -> None:
        super().__init__()

    @abstractclassmethod
    def get_matches(cls, ligand:Ligand) -> Substructure:
        pass