from rdkit import Chem
from rdkit.Chem import Mol
from AaronTools.fileIO import FileReader
from .pattern import Pattern
from AaronTools.geometry import Geometry
from enum import Enum, IntEnum

def test_function():
    print('This is a test function nice')


class Coordinates:
    def __init__(self, coordinates:tuple):
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]


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

class Atom:
    def __init__(self, element:Element, bonds:list['Bond']=None, coordinates:Coordinates=None):
        self.element:Element = element

        self.bonds:list['Bond'] = []
        for bond in bonds:
            if isinstance(bond, 'Bond'):
                self.bonds.append(bond)
            else:
                raise ValueError
        
        self.coordinates:Coordinates = None
        if isinstance(coordinates, Coordinates):
            self.coodinates = coordinates
        else:
            raise ValueError()

class Bond:
    def __init__(self, atom1, atom2, multiplicity:int=1):
        self.atom1:Atom = atom1
        self.atom2:Atom = atom2
        self.multiplicity:int = multiplicity

    @property
    def atoms(self) -> tuple:
        return (self.atom1, self.atom2)



class Denticity(IntEnum):
    Monodentate  = 1
    Bidentate    = 2
    Tridentate   = 3
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