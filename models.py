# models.py
import math

class Location:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def distance_to(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return (dx * dx + dy * dy + dz * dz) ** 0.5

class Atom:
    def __init__(self, atom_id, name, atom_type, charge, location, temp_factor, mer=None):
        self.atom_id = atom_id
        self.name = name
        self.type = atom_type
        self.charge = charge
        self.location = location
        self.temp_factor = temp_factor
        self.mer = mer

class Mer:
    def __init__(self, name, mer_type, mer_id, chain):
        self.name = name
        self.type = mer_type
        self.mer_id = mer_id
        self.chain = chain
        self.atoms = []
        self.bond_count = {}  # Mer -> int
        self.center_of_mass = Location(0,0,0)

    def add_atom(self, atom):
        self.atoms.append(atom)
        atom.mer = self

    def calc_com(self):
        if not self.atoms:
            self.center_of_mass = Location(0,0,0)
            return
        x = sum(a.location.x for a in self.atoms)
        y = sum(a.location.y for a in self.atoms)
        z = sum(a.location.z for a in self.atoms)
        count = len(self.atoms)
        self.center_of_mass = Location(x/count, y/count, z/count)

    def add_bond(self, other):
        self.bond_count[other] = self.bond_count.get(other, 0) + 1

    def is_bonded_to(self, other, threshold=4.0):
        # Check distance between any atom pairs
        bonded = False
        for a1 in self.atoms:
            for a2 in other.atoms:
                dist = math.sqrt((a1.location.x - a2.location.x)**2 + 
                                 (a1.location.y - a2.location.y)**2 +
                                 (a1.location.z - a2.location.z)**2)
                if dist <= threshold:
                    self.add_bond(other)
                    other.add_bond(self)
                    bonded = True
                    break
            if bonded:
                break
        return bonded

    def get_bond_count_with(self, other):
        return self.bond_count.get(other, 0)

    def __eq__(self, other):
        if not isinstance(other, Mer):
            return False
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


class Interaction:
    def __init__(self, from_mer, to_mer):
        self.from_mer = from_mer.name
        self.to_mer = to_mer.name
        bond_count = from_mer.get_bond_count_with(to_mer)
        # Avoid division by zero
        from_size = len(from_mer.atoms) or 1
        to_size = len(to_mer.atoms) or 1
        self.affinity = bond_count / math.sqrt(from_size * to_size)

    def __str__(self):
        return f"{self.from_mer} --[{self.affinity}]--> {self.to_mer}"
