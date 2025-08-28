import os
import json
from abc import ABC, abstractmethod
from .ligand import Element, Coordinate, Atom, Bond, Ligand, MonodentateLigand, BidentateLigand, TridentateLigand, TetradentateLigand



class LigandParser(ABC):
    
    def __init__(self) -> None:
        self._parsed_ligands:dict[str, Ligand] = {}
        self._filepaths:dict[str, str] = {}

    @staticmethod
    def get_ccdc_code(filepath) -> str:
        return os.path.basename(filepath)[0:6]

    #Parse either an sdf or mol file
    #Cannot just use aarontools parser since doesn't parse comment lines. But can rip code
    def parse_legacy(self, filepath:str) -> Ligand:
        name = self.get_ccdc_code(filepath)

        ligand = None
        atoms:dict[int,Atom] = {}
        bonds:list[Bond] = []
        charge = None
        coordinating_atoms:list[Atom] = None
        formal_charge_dict = {}

        lines = []
        with open(filepath, 'r') as inf:
            lines = inf.readlines()

        for i, line in enumerate(lines):
            if i == 0:
                header = line.strip()
                coordinating_atom_indexes = tuple(header.split(','))
                denticity = len(coordinating_atom_indexes)

            if i == 2:
                comment = line.strip()
                charge = int(comment.split()[:1])

            if i == 3:
                natoms = int(line[0:3])
                nbonds = int(line[3:6])

            if i == 4:
                atoms = {}
                atom_index = 1
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
                    valence          = int(line[48:51].strip())

                    
                    element = Element[element_str]
                    coordinate = Coordinate((x_coord, y_coord, z_coord))
                    atom = Atom()
                    atom.set_element(element=element)
                    atom.set_coordinate(coordinate=coordinate)
                    atom.set_valence(valence=valence)

                    atoms[atom_index] = atom
                    atom_index += 1
                    
                bonds = []
                for line in lines[i + natoms : i + natoms + nbonds]:
                    a1_index = int(line[0:3].strip())
                    a2_index = int(line[3:6].strip())
                    m  = int(line[6:9].strip())

                    a1 = atoms[a1_index]
                    a2 = atoms[a2_index]
                    bond = Bond(a1, a2, m)
                    bonds.append(bond)
                    a1.add_bond(bond)
                    a2.add_bond(bond)

                match (denticity):
                    case 1:
                        ligand = MonodentateLigand()
                    case 2:
                        ligand = BidentateLigand()
                    case 3:
                        ligand = TridentateLigand()
                    case 4:
                        ligand = TetradentateLigand()
                
                for line in lines[i + natoms + nbonds:]:
                    if "M  CHG" in line:
                        formal_charge_line = line
                        formal_charge_list = [int(i) for i in formal_charge_line.split()[3:]]
                        for i in range(0,len(formal_charge_list),2):
                            formal_charge_dict[formal_charge_list[i]] = formal_charge_list[i+1]
        
        for index, charge in formal_charge_dict.items():
            atoms[index].set_formal_charge(charge)

        ligand.set_name(name=name)
        ligand.set_charge(charge=charge)
        for atom in atoms:
            ligand.add_atom(atom)
        for bond in bonds:
            ligand.add_bond(bond)

        coordinating_atoms = [atoms[i] for i in coordinating_atom_indexes]
        ligand.set_coordinating_atoms(coordinating_atoms=coordinating_atoms)
        
        self._filepaths[ligand.name] = filepath
        self._parsed_ligands[ligand.name] = ligand


    def parse(self, filepath):
        raise NotImplementedError
        
        lines = []
        with open(filepath, 'r') as inf:
            lines = inf.readlines()
        
        for i, line in enumerate(lines):
            pass

        
    def __getitem__(self, key):
        return self._parsed_ligands[key]

    def get_filepath(self, ligand: str|Ligand):
        if isinstance(ligand, Ligand):
            return self._filepaths[ligand.name]
        else:
            return self._filepaths[ligand]
    
class LigandWriter:
    def __init__(self):
        pass

    #TODO: refactor so not duplicating so much code.
    @staticmethod
    def write_mol(ligand:Ligand, filepath):
        with open(filepath, 'w') as outf:
            outf.write(ligand.name + '\n')
            outf.write('Generated By RL-Tools\n')
            header_data = {
                "name" : ligand.name,
                "charge" : ligand.charge,
                "connections" : ligand.get_connection_indexes()
            }
            json_string = json.dumps(header_data)
            outf.write(json_string + '\n')
            
            atom_count = len(ligand.atoms)
            bond_count = len(ligand.bonds)
            counts_line = f'{atom_count.rjust(3)}{bond_count.rjust(3)}  0  0  0  0  0  0  0  0999 V2000'
            outf.write(counts_line + '\n')
            for atom in ligand.atoms:
                atom_line = atom.get_file_str()
                outf.write(atom_line + '\n')
            
            for bond in ligand.bonds:
                bond_line = bond.get_file_str()
                outf.write(bond_line + '\n')

            formal_charge_line = ligand.get_formal_charge_line()
            if formal_charge_line:
                outf.write(formal_charge_line + '\n')
            
            outf.write('M  END')
    
    @staticmethod
    def write_sdf(ligands:Ligand|list[Ligand], filepath, append=True):
        if isinstance(ligands, Ligand):
            LigandWriter.write_mol(ligands, filepath)
            return
        
        write_type = 'a' if append == True else 'w'
        with open(filepath, write_type) as outf:
            for ligand in ligands:
                outf.write(ligand.name + '\n')
                outf.write('Generated By RL-Tools\n')
                header_data = {
                    "name" : ligand.name,
                    "charge" : ligand.charge,
                    "connections" : ligand.get_connection_indexes()
                }
                json_string = json.dumps(header_data)
                outf.write(json_string + '\n')
                
                atom_count = len(ligand.atoms)
                bond_count = len(ligand.bonds)
                counts_line = f'{atom_count.rjust(3)}{bond_count.rjust(3)}  0  0  0  0  0  0  0  0999 V2000'
                outf.write(counts_line + '\n')
                for atom in ligand.atoms:
                    atom_line = atom.get_file_str()
                    outf.write(atom_line + '\n')
                
                for bond in ligand.bonds:
                    bond_line = bond.get_file_str()
                    outf.write(bond_line + '\n')

                formal_charge_line = ligand.get_formal_charge_line()
                if formal_charge_line:
                    outf.write(formal_charge_line + '\n')
                
                outf.write('M  END\n')
                outf.write('\n')
                outf.write('$$$$\n')
