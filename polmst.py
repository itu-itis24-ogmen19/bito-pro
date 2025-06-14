import math
import logging
import heapq
from io import StringIO
from models import Atom, Mer, Location, Interaction

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# ----------------------------------------------------------------------
# Parsing helpers (unchanged)
# ----------------------------------------------------------------------
def read_input_file(file_path):
    with open(file_path, 'r') as f:
        return f.readlines()

def parse_atoms(raw_lines):
    mers = {}
    for line in raw_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = parse_atom_line(line)
            if atom:
                mer_name = atom.mer.name
                if mer_name not in mers:
                    mers[mer_name] = atom.mer
                mers[mer_name].add_atom(atom)

    for mer in mers.values():
        mer.calc_com()

    return mers

def parse_atom_line(line):
    try:
        if len(line) < 66:
            return None

        atom_id = int(line[6:11].strip())
        name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_number = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        temp_factor = float(line[60:66].strip())

        location = Location(x, y, z)
        mer_name = f"{residue_name}-{residue_number}({chain_id})"
        mer = Mer(mer_name, residue_name, residue_number, chain_id)
        atom = Atom(atom_id, name, residue_name, "", location, temp_factor, mer)
        return atom
    except:
        return None

# ----------------------------------------------------------------------
# Graph construction helpers (unchanged)
# ----------------------------------------------------------------------
def calculate_interactions(mers):
    mer_list = list(mers.values())

    for i in range(len(mer_list)):
        src = mer_list[i]
        for j in range(len(mer_list)):
            if i == j:
                continue
            dst = mer_list[j]
            bonded = False
            for src_atom in src.atoms:
                for dst_atom in dst.atoms:
                    if src_atom.location.distance_to(dst_atom.location) < 4.5:
                        src.add_bond(dst)
                        dst.add_bond(src)
                        bonded = True
                        break
                if bonded:
                    break

    interactions = []
    for i in range(len(mer_list)):
        for j in range(i + 1, len(mer_list)):
            src = mer_list[i]
            dst = mer_list[j]
            if src.is_bonded_to(dst):
                inter = Interaction(src, dst)
                interactions.append(inter)

    return interactions

def build_adjacency_map(mers):
    adjacency_map = {m.name: {} for m in mers.values()}

    for from_mer in mers.values():
        for to_mer_name, bond_count in from_mer.bond_count.items():
            if bond_count > 0 and to_mer_name in mers:
                to_mer = mers[to_mer_name]
                affinity = bond_count / math.sqrt(len(from_mer.atoms) * len(to_mer.atoms))
                weight = 1.0 / affinity
                adjacency_map[from_mer.name][to_mer_name] = weight
                adjacency_map[to_mer_name][from_mer.name] = weight

    return adjacency_map

def dijkstra(adjacency_map, source_mer):
    distances = {m: math.inf for m in adjacency_map}
    distances[source_mer] = 0.0
    visited = set()
    unvisited = set(adjacency_map.keys())

    while unvisited:
        current = min(unvisited, key=lambda m: distances[m])
        unvisited.remove(current)
        visited.add(current)

        if distances[current] == math.inf:
            break

        for neighbor, weight in adjacency_map[current].items():
            if neighbor not in visited:
                new_dist = distances[current] + weight
                if new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist

    return distances

def find_connected_components(adjacency_map):
    visited = set()
    components = []

    for mer_name in adjacency_map:
        if mer_name not in visited:
            stack = [mer_name]
            comp = set()
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    comp.add(node)
                    for nbr in adjacency_map[node]:
                        if nbr not in visited:
                            stack.append(nbr)
            components.append(comp)

    return components

def prune_to_component(adjacency_map, keep_set):
    new_map = {}
    for m in keep_set:
        if m in adjacency_map:
            sub_nbrs = {nbr: w for nbr, w in adjacency_map[m].items() if nbr in keep_set}
            new_map[m] = sub_nbrs
    return new_map

def prune_interactions(interactions, sub_map):
    return [i for i in interactions if i.from_mer in sub_map and i.to_mer in sub_map]

def prune_mers(mers, sub_map):
    return {m: mers[m] for m in sub_map if m in mers}

# ----------------------------------------------------------------------
# *** UPDATED FUNCTION ***
# ----------------------------------------------------------------------
def generate_enhanced_pdb(weights, source_mer, original_file, mers):
    """
    Return a PDB that keeps every original atom record so the
    downloaded file’s atom count matches the input exactly.

    For each atom we overwrite its B-factor (TempFactor column)
    with the weight/dist-to-source of its parent Mer, allowing
    external tools to colour by that value.

    A single REMARK (type 999) is added just before the first
    coordinate line to document what the B-factor encodes.
    """
    out = StringIO()

    with open(original_file, "r") as pdb_in:
        wrote_remark = False

        for line in pdb_in:
            is_atom = line.startswith("ATOM") or line.startswith("HETATM")

            # ----------------------------------------------------------
            # Insert one explanatory remark before the first ATOM line
            # ----------------------------------------------------------
            if is_atom and not wrote_remark:
                out.write(f"REMARK 999  Distances from {source_mer} are stored in the B-factor column\n")
                wrote_remark = True

            # ----------------------------------------------------------
            # If this is a coordinate line, patch columns 61-66
            # ----------------------------------------------------------
            if is_atom and len(line) >= 66:
                resname = line[17:20].strip()
                chain   = line[21].strip()
                resnum  = line[22:26].strip()
                mer_name = f"{resname}-{resnum}({chain})"

                new_b = weights.get(mer_name, None)
                if new_b is None or math.isinf(new_b) or math.isnan(new_b):
                    # keep original record untouched
                    out.write(line)
                else:
                    new_b = min(new_b, 999.99)  # fit into %6.2f
                    patched = f"{line[:60]}{new_b:6.2f}{line[66:]}"
                    out.write(patched)
            else:
                out.write(line)

    return out.getvalue()

# ----------------------------------------------------------------------
# Main processing pipeline (unchanged)
# ----------------------------------------------------------------------
def process_pdb_file(file_path):
    raw_lines = read_input_file(file_path)
    mers = parse_atoms(raw_lines)

    interactions = calculate_interactions(mers)
    adjacency_map = build_adjacency_map(mers)

    comps = find_connected_components(adjacency_map)
    largest_comp = max(comps, key=len)

    adjacency_map = prune_to_component(adjacency_map, largest_comp)
    interactions  = prune_interactions(interactions, adjacency_map)
    mers          = prune_mers(mers, adjacency_map)

    best_source_mer = None
    min_sum_dist = math.inf
    all_distance_sums = {}

    for mer_name in largest_comp:
        dist_map = dijkstra(adjacency_map, mer_name)
        total_dist = sum(d for d in dist_map.values() if not math.isinf(d))
        all_distance_sums[mer_name] = dist_map
        if total_dist < min_sum_dist:
            min_sum_dist = total_dist
            best_source_mer = mer_name

    distances_best = all_distance_sums[best_source_mer]
    total_weight_sums = dict(distances_best)

    return (
        best_source_mer,
        mers,
        total_weight_sums,
        interactions,
        all_distance_sums,
    )
