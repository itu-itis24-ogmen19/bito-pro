import math

class Location:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def distance_to(self, other):
        return math.sqrt(
            (self.x - other.x) ** 2 +
            (self.y - other.y) ** 2 +
            (self.z - other.z) ** 2
        )

class Atom:
    def __init__(self, atom_id, name, residue_name, element, location, temp_factor, mer):
        self.atom_id = atom_id
        self.name = name
        self.residue_name = residue_name
        self.element = element
        self.location = location
        self.temp_factor = temp_factor
        self.mer = mer  # Mer object

class Mer:
    def __init__(self, name, residue_name, residue_number, chain):
        self.name = name
        self.residue_name = residue_name
        self.residue_number = residue_number
        self.chain = chain
        self.atoms = []
        self.bond_count = {}  # dict of mer_name: bond_count
        self.center_of_mass = None

    def add_atom(self, atom):
        self.atoms.append(atom)

    def calc_com(self):
        if not self.atoms:
            self.center_of_mass = None
            return
        x = sum(atom.location.x for atom in self.atoms) / len(self.atoms)
        y = sum(atom.location.y for atom in self.atoms) / len(self.atoms)
        z = sum(atom.location.z for atom in self.atoms) / len(self.atoms)
        self.center_of_mass = Location(x, y, z)

    def add_bond(self, other_mer):
        if other_mer.name in self.bond_count:
            self.bond_count[other_mer.name] += 1
        else:
            self.bond_count[other_mer.name] = 1

    def is_bonded_to(self, other_mer):
        return other_mer.name in self.bond_count and self.bond_count[other_mer.name] > 0

    def get_bond_count_with(self, other_mer):
        return self.bond_count.get(other_mer.name, 0)

class Interaction:
    def __init__(self, from_mer, to_mer):
        """
        Interaction between two Mers, storing their names and the calculated weight.
        Weight is the primary metric used for all calculations and visualizations.
        """
        self.from_mer = from_mer.name  # string
        self.to_mer = to_mer.name      # string
        bond_count = from_mer.get_bond_count_with(to_mer)
        # Avoid division by zero
        from_size = len(from_mer.atoms) or 1
        to_size = len(to_mer.atoms) or 1
        if bond_count > 0:
            self.weight = 1.0 / (bond_count / math.sqrt(from_size * to_size))  # weight = 1 / affinity
        else:
            self.weight = math.inf  # Represents no interaction

    def __str__(self):
        return f"{self.from_mer} --[{self.weight:.2f}]--> {self.to_mer}"
