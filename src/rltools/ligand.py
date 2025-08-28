import os
from rdkit import Chem
from rdkit.Chem import Mol
from AaronTools.fileIO import FileReader
from .pattern import Pattern
from AaronTools.geometry import Geometry
from enum import Enum, IntEnum
from abc import ABC, abstractmethod
import json

class Element(IntEnum):
    # --- Real elements ---
    H  = 1
    He = 2
    Li = 3
    Be = 4
    B  = 5
    C  = 6
    N  = 7
    O  = 8
    F  = 9
    Ne = 10
    Na = 11
    Mg = 12
    Al = 13
    Si = 14
    P  = 15
    S  = 16
    Cl = 17
    Ar = 18
    K  = 19
    Ca = 20
    Sc = 21
    Ti = 22
    V  = 23
    Cr = 24
    Mn = 25
    Fe = 26
    Co = 27
    Ni = 28
    Cu = 29
    Zn = 30
    Ga = 31
    Ge = 32
    As = 33
    Se = 34
    Br = 35
    Kr = 36
    Rb = 37
    Sr = 38
    Y  = 39
    Zr = 40
    Nb = 41
    Mo = 42
    Tc = 43
    Ru = 44
    Rh = 45
    Pd = 46
    Ag = 47
    Cd = 48
    In = 49
    Sn = 50
    Sb = 51
    Te = 52
    I  = 53
    Xe = 54
    Cs = 55
    Ba = 56
    La = 57
    Ce = 58
    Pr = 59
    Nd = 60
    Pm = 61
    Sm = 62
    Eu = 63
    Gd = 64
    Tb = 65
    Dy = 66
    Ho = 67
    Er = 68
    Tm = 69
    Yb = 70
    Lu = 71
    Hf = 72
    Ta = 73
    W  = 74
    Re = 75
    Os = 76
    Ir = 77
    Pt = 78
    Au = 79
    Hg = 80
    Tl = 81
    Pb = 82
    Bi = 83
    Po = 84
    At = 85
    Rn = 86
    Fr = 87
    Ra = 88
    Ac = 89
    Th = 90
    Pa = 91
    U  = 92
    Np = 93
    Pu = 94
    Am = 95
    Cm = 96
    Bk = 97
    Cf = 98
    Es = 99
    Fm = 100
    Md = 101
    No = 102
    Lr = 103
    Rf = 104
    Db = 105
    Sg = 106
    Bh = 107
    Hs = 108
    Mt = 109
    Ds = 110
    Rg = 111
    Cn = 112
    Nh = 113
    Fl = 114
    Mc = 115
    Lv = 116
    Ts = 117
    Og = 118

    # --- Pseudo-elements / placeholders ---
    Any =  0      # "*" wildcard atom
    Du  = -1     # Dummy atom
    R   = -2     # Generic R group
    R1  = -3
    R2  = -4
    R3  = -5
    R4  = -6
    X   = -7     # Generic halogen

class Coordinate:
    def __init__(self, coordinates:tuple):
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]
    
    def __hash__(self):
        return hash((self.x, self.y, self.z))
    
    def __eq__(self, value):
        if isinstance(value, Coordinate):
            return self.x == value.x and self.y == value.y and self.z == value.z
        return False

class Atom:
    def __init__(self):
        self.element:Element = None
        self.coordinate:Coordinate = None
        self.bonds:set[Bond] = set()
        self.valence = None
        self.formal_charge = None

    def set_element(self, element:Element):
        if isinstance(element, Element):
            self.element = element
        else:
            raise ValueError()

    def set_coordinate(self, coordinate:Coordinate):
        if isinstance(coordinate, Coordinate):
            self.coordinate = Coordinate
        else:
            raise ValueError()
        
    def set_valence(self, valence:int):
        self.valence = valence

    def set_formal_charge(self, charge):
        self.formal_charge = charge

    def add_bond(self, bond:'Bond'):
        if not isinstance(bond, Bond):
            raise ValueError()
        
        if self not in bond.atoms:
            raise ValueError()
        
        self.bonds.add(bond)
        
    def get_file_str(self):
        x = str(round(self.coordinate.x, 5))
        y = str(round(self.coordinate.y, 5))
        z = str(round(self.coordinate.z, 5))
        element_str = self.element.name + (" " if len(self.element.name) == 1 else "")
        valence_str = str(self.valence)

        return f"{x.rjust(10)}{y.rjust(10)}{z.rjust(10)}{element_str.rjust(3)}  0  0  0  0  0{valence_str.rjust(3)}  0  0  0  0  0  0"

    def __hash__(self):
        return hash((self.coordinate, self.element))
    
    def __eq__(self, value):
        if not isinstance(value, Atom):
            return False
        if self.coordinate != value.coordinate:
            return False
        
        # Including this results in infinite recursion
        # if self.bonds != value.bonds:
        #     return False

        return True
        
    
    
    

class Bond:
    def __init__(self, atom1, atom2, multiplicity:int=1):
        self.atom1:Atom = atom1
        self.atom2:Atom = atom2
        self.multiplicity:int = multiplicity
    
    def __hash__(self):
        return hash((frozenset((self.atom1, self.atom2)), self.multiplicity))
    
    def __eq__(self, value):
        if not isinstance(value, Bond):
            return False
        
        if frozenset(self.atoms) != frozenset(value.atoms):
            return False

        return True

    @property
    def atoms(self) -> tuple:
        return (self.atom1, self.atom2)

    def get_file_str(self, atom_list:list[Atom]):
        a1_position = str(atom_list.index(self.atom1))
        a2_position = str(atom_list.index(self.atom2))
        m = str(self.multiplicity)
        return f'{a1_position.rjust(3)}{a2_position.rjust(3)}{m.rjust(3)}'



class Denticity(IntEnum):
    Monodentate  = 1
    Bidentate    = 2
    Tridentate   = 3
    Tetradentate = 4

class Ligand(ABC):
    denticity:Denticity = None

    def __init__(self) -> None:
        super().__init__()
        self.name = ""
        self.charge:int = None
        self.atoms:list[Atom] = []
        self.bonds:list[Bond] = []
        self.coordinating_atoms:list[Atom] = []

    def __hash__(self):
        return hash((frozenset(self.atoms)))
    
    def __eq__(self, value: object) -> bool:
        if not isinstance(value, self.__class__):
            return False
        if self.charge != value.charge:
            return False
        if frozenset(self.atoms) != frozenset(value.atoms):
            return False
        # if frozenset(self.bonds) != frozenset(value.bonds):
        #     return False
        # if frozenset(self.coordinating_atoms) != frozenset(self.coordinating_atoms):
        #     return False

        return True


    def set_name(self, name:str) -> None:
        self.name = name
    
    def set_charge(self, charge:int) -> None:
        self.charge = charge
    
    def set_coordinating_atoms(self, coordinating_atoms:list[Atom]):
        for atom in coordinating_atoms:
            if atom not in self.atoms:
                raise ValueError('atom cannot be coordinating atom if not a part of the ligand')
        
        self.coordinating_atoms = coordinating_atoms

    def add_atom(self, atom:Atom) -> list[Atom]:
        self.atoms.append(atom)
        return self.atoms

    def add_bond(self, bond:Bond) -> list[Bond]:
        self.bonds.append(bond)
        return self.bonds

    def match_pattern(self, pattern:Pattern):
        matches = pattern.get_matches(self)
        return matches
    



class MonodentateLigand(Ligand):
    denticity = Denticity(1)

class BidentateLigand(Ligand):
    denticity = Denticity(2)

class TridentateLigand(Ligand):
    denticity = Denticity(3)
    
    def __init__(self):
        super().__init__()

    @property
    def side_connections(self) -> tuple[Atom]:
        pass
        #logic comparing atom distances
    
class TetradentateLigand(Ligand):
    denticity = Denticity(4)




                
                


            






if __name__ == "__main__":
    ligand = Ligand()