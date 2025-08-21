from rdkit import Chem
from rdkit.Chem import Mol
from AaronTools.fileIO import FileReader
from .pattern import Pattern
from AaronTools.geometry import Geometry
from enum import Enum, IntEnum
from abc import ABC, abstractmethod

def test_function():
    print('This is a test function nice')


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
    def __init__(self):
        self.element:Element = None
        self.coordinate:Coordinate = None
        self.bonds:set[Bond] = set()

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

    def add_bond(self, bond:'Bond'):
        if not isinstance(bond, 'Bond'):
            raise ValueError()
        
        if self not in bond.atoms:
            raise ValueError()
        
        self.bonds.add(bond)
        
    def __hash__(self):
        return hash((self.coordinate, self.element))
    
    def __eq__(self, value):
        if not isinstance(value, Atom):
            return False
        if self.coordinate != value.coordinate:
            return False
        if self.bonds != value.bonds:
            return False

        return True
        
    
    
    

class Bond:
    def __init__(self, atom1, atom2, multiplicity:int=1):
        self.atom1:Atom = atom1
        self.atom2:Atom = atom2
        self.multiplicity:int = multiplicity

    @property
    def atoms(self) -> tuple:
        return (self.atom1, self.atom2)
    
    def __hash__(self):
        return hash((frozenset(self.atom1, self.atom2), self.multiplicity))
    
    def __eq__(self, value):
        if not isinstance(value, Bond):
            return False
        
        if frozenset(self.atoms) != frozenset(value.atoms):
            return False

        return True
        



class Denticity(IntEnum):
    Monodentate  = 1
    Bidentate    = 2
    Tridentate   = 3
    Tetradentate = 4

class Ligand(ABC):
    denticity:Denticity = None

    def __init__(self) -> None:
        super.__init__()
        self.name = ""
        self.charge:int = None
        self.atoms:list[Atom] = []
        self.bonds:list[Bond] = []
        self.coordinating_atoms:list[Atom] = []


    def set_name(self, name:str) -> None:
        self.name = name
    
    def set_charge(self, charge:int) -> None:
        self.charge = charge
    
    def set_coordinating_atoms(self, coordinating_atoms:list[Atom]):
        for atom in coordinating_atoms:
            if atom not in self.atoms:
                raise ValueError('atom cannot be coordinating atom if not a part of the ligand')
        
        self.coordinating_atoms = coordinating_atoms

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


class LigandParser(ABC):
    
    def __init__(self) -> None:
        self._parsed_ligands:dict[str, Ligand] = {}
        self._filepaths:dict[str, str] = {}


    #Parse either an sdf or mol file
    #Cannot just use aarontools parser since doesn't parse comment lines. But can rip code
    def parse(self, filepath:str) -> Ligand:
        header = None
        natoms = None
        nbonds = None
        atoms = []

        with open(filepath, 'r') as inf:
            lines = inf.readlines()
            self.all_geom = []
            progress = 0
            for i, line in enumerate(lines):
                progress += 1
                if "$$$$" in line:
                    progress = 0
                    match (denticity):
                        case 1:
                            ligand = MonodentateLigand()
                        case 2:
                            ligand = BidentateLigand()
                        case 3:
                            ligand = TridentateLigand()
                        case 4:
                            ligand = TetradentateLigand()

                    ligand.set_name(name=name)
                    ligand.set_charge(charge=charge)

                    coordinating_atoms = [atoms[i] for i in coordinating_atom_indexes]
                    ligand.set_coordinating_atoms(coordinating_atoms=coordinating_atoms)

                    continue
                
                if progress == 1:
                    header = line.strip()
                    coordinating_atom_indexes = tuple(header.split(','))

                if progress == 3:
                    comment = line.strip()
                    charge = int(comment.split()[:1])

                if progress == 4:
                    natoms = int(line[0:3])
                    nbonds = int(line[3:6])

                if progress == 5:
                    atoms = {}
                    index = 1
                    for line in lines[i : i + natoms]:
                        
                        x_coord          = int(line[0:10].strip())
                        y_coord          = int(line[10:20].strip())
                        z_coord          = int(line[20:30].strip())
                        element_str      =     line[30:34].strip()
                        mass_diff        = int(line[34:36].strip())
                        charge_code      = int(line[36:39].strip())
                        stereo_parity    = int(line[39:42].strip())
                        h0_designator    = int(line[42:45].strip())
                        stereo_care_flag = int(line[45:48].strip())
                        
                        element = Element[element_str]
                        coordinate = Coordinate((x_coord, y_coord, z_coord))
                        atom = Atom()
                        atom.set_element(element=element)
                        atom.set_coordinate(coordinate=coordinate)

                        atoms[index] = atom
                        index += 1
                        
                    bonds = []
                    try:
                        for line in lines[i + natoms : i + natoms + nbonds]:
                            a1 = int(line[0:3].strip())
                            a2 = int(line[3:6].strip())
                            m  = int(line[])
                    except ValueError:
                        for line in lines[i + natoms : i + natoms + nbonds]:
                            a1, a2, _ = line.split()
                            a1 = int(a1)
                            a2 = int(a2)
                            self.atoms[a1].connected.add(self.atoms[a2])
                            self.atoms[a2].connected.add(self.atoms[a1])

                    for j, a in enumerate(self.atoms):
                        a.name = str(j + 1)

                    self.other["charge"] = 0
                    for line in lines[i + natoms + nbonds:]:
                        if "CHG" in line:
                            self.other["charge"] += int(line.split()[-1])
                        if "$$$$" in line:
                            break


    def __getitem__(self, key):
        return self._parsed_ligands[key]

    def get_filepath(self, ligand):
        if isinstance(ligand, Ligand):
            return self._filepaths[ligand.name]
        else:
            return self._filepaths[ligand]
    


if __name__ == "__main__":
    ligand = Ligand()