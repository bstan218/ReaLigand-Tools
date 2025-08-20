from rdkit import Chem
from rdkit.Chem import Mol
from AaronTools.fileIO import FileReader
from .pattern import Pattern
from AaronTools.geometry import Geometry
from enum import Enum

def test_function():
    print('This is a test function nice')


class Bond:
    pass

class Atom:
    pass

class Denticity(Enum):
    Monodentate = 1
    Bidentate = 2
    Tridentate = 3
    Tetradentate = 4

class Ligand:
    def __init__(self, filepath:str) -> None:
        self.charge:int = None
        self.denticity:Denticity = None
        self.atoms:list[Atom] = None
        self.bonds:list[Bond] = None
        
        
        filereader = FileReader(filepath)
        
        self.rdkit_mol:Mol = Chem.MolFromMolFile(filepath, removeHs=False)
        output = filereader.read_file(get_all=True)
        print(output)

    def match_pattern(self, pattern:Pattern):
        matches = pattern.get_matches(self)
        return matches


class LigandParser:
    def __init__(self) -> None:
        pass

    def parse(self, filepath:str) -> Ligand:
        rdkit_mol = Chem.MolFromMolFile(filepath, removeHs=False)
        aaron_filereader = FileReader(filepath)
        aaron_geom = aaron_filereader.read_file()
        ligand = Ligand()

if __name__ == "__main__":
    ligand = Ligand()