import math
import logging
from models import Atom, Mer, Location, Interaction

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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

def calculate_interactions(mers):
    interactions = []
    mer_list = list(mers.values())
    for i in range(len(mer_list)):
        for j in range(i+1, len(mer_list)):
            src = mer_list[i]
            dst = mer_list[j]
            if src.is_bonded_to(dst):
                inter = Interaction(src, dst)
                interactions.append(inter)
    return interactions

def build_adjacency_map(mers):
    adjacency_map = {}
    for mer in mers.values():
        adjacency_map[mer] = {}

    all_mers = list(mers.values())
    for i, from_mer in enumerate(all_mers):
        for j, to_mer in enumerate(all_mers):
            if i == j:
                continue
            bond_count = from_mer.get_bond_count_with(to_mer)
            if bond_count > 0:
                affinity = bond_count / math.sqrt(len(from_mer.atoms) * len(to_mer.atoms))
                weight = 1.0 / affinity
                adjacency_map[from_mer][to_mer] = weight
    return adjacency_map

def dijkstra(adjacency_map, source_mer):
    distances = {m: math.inf for m in adjacency_map.keys()}
    distances[source_mer] = 0.0
    visited = set()
    unvisited = list(adjacency_map.keys())

    while unvisited:
        current = min(unvisited, key=lambda m: distances[m])
        unvisited.remove(current)
        visited.add(current)

        for neighbor, weight in adjacency_map[current].items():
            if neighbor in visited:
                continue
            new_dist = distances[current] + weight
            if new_dist < distances[neighbor]:
                distances[neighbor] = new_dist
    return distances

def find_best_source_mer(mers, adjacency_map):
    best_mer = None
    min_total_distance = math.inf
    for source_mer in mers:
        dist = dijkstra(adjacency_map, source_mer)
        total_distance = sum(d for d in dist.values() if not math.isinf(d))
        if total_distance < min_total_distance:
            min_total_distance = total_distance
            best_mer = source_mer
    return best_mer

def format_atom_line(atom):
    recordName = "ATOM  "
    serial = atom.atom_id
    name = f"{atom.name:<4}"
    resName = f"{atom.type:<3}"
    chainID = atom.mer.chain
    resSeq = atom.mer.mer_id
    x = atom.location.x
    y = atom.location.y
    z = atom.location.z
    occupancy = 1.0
    tempFactor = atom.temp_factor

    line = f"{recordName:6s}{serial:5d} {name:>4s} {resName:>3s} {chainID:1s}{resSeq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempFactor:6.2f}"
    return line

def generate_enhanced_pdb(raw_lines, source_mer, shortest_paths, interactions, mers):
    # Compute maxDistance
    finite_distances = [d for d in shortest_paths.values() if not math.isinf(d)]

    if not finite_distances:
        # No finite distances, default to 1.0
        max_distance = 1.0
    else:
        max_distance = max(finite_distances)
        # If max_distance is 0.0, set it to 1.0 to prevent division by zero
        if max_distance == 0.0:
            max_distance = 1.0

    # Update tempFactor of atoms based on distances
    for mer in mers.values():
        distance = shortest_paths.get(mer, math.inf)
        temp_factor_value = 0.0 if math.isinf(distance) else (distance / max_distance) * 100.0
        for atom in mer.atoms:
            atom.temp_factor = temp_factor_value

    sb = []
    # Append non-ATOM/HETATM lines
    for line in raw_lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            sb.append(line.strip())

    # Append updated ATOM lines
    for mer in mers.values():
        for atom in mer.atoms:
            atom_line = format_atom_line(atom)
            sb.append(atom_line)

    sb.append(f"REMARK 1 Source Mer: {source_mer.name}")
    sb.append("REMARK 2 Shortest Paths:")
    for m, d in shortest_paths.items():
        sb.append(f"REMARK 3 {m.name}: {d}")
    sb.append("REMARK 4 Interactions:")
    for inter in interactions:
        sb.append(f"REMARK 5 {inter}")
    return "\n".join(sb)

def process_pdb_file(file_path):
    raw_lines = read_input_file(file_path)
    mers = parse_atoms(raw_lines)
    interactions = calculate_interactions(mers)
    adjacency_map = build_adjacency_map(mers)
    best_source_mer = find_best_source_mer(mers.values(), adjacency_map)

    # Generate enhanced PDB for every mer as the source, not just best source mer
    all_enhanced_pdbs = {}
    all_distances = {}
    for source_mer in mers.values():
        distances = dijkstra(adjacency_map, source_mer)
        enhanced_pdb = generate_enhanced_pdb(raw_lines, source_mer, distances, interactions, mers)
        all_enhanced_pdbs[source_mer.name] = enhanced_pdb
        all_distances[source_mer.name] = distances

    # shortest_paths for the best_source_mer:
    shortest_paths = all_distances[best_source_mer.name]

    # Now we return all enhanced pdb files so that the UI can provide download options for each mer
    return best_source_mer.name, mers, shortest_paths, all_enhanced_pdbs
