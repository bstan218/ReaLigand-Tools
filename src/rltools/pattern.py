from __future__ import annotations
from typing import TYPE_CHECKING
from abc import ABC, abstractclassmethod

if TYPE_CHECKING:
    from .ligand import Ligand

class Substructure:
    def __init__(self) -> None:
        pass

class Pattern(ABC):
    def __init__(self) -> None:
        super().__init__()

    @abstractclassmethod
    def get_matches(cls, ligand:Ligand) -> Substructure:
        pass



class Cycle(Pattern):
    def get_matches(cls, ligand:Ligand) -> Substructure:
        pass

class Aromatic(Pattern):
    pass

class Bridge(Pattern):
    pass


class RGroup(Pattern):
    pass

class ExtendedRGroup(Pattern):
    pass

class ArylExtendedRGroup(Pattern):
    pass

class AlkylExtendedRGroup(Pattern):
    pass