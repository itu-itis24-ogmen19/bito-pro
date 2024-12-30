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
    adjacency_map = {}
    for mer in mers.values():
        adjacency_map[mer.name] = {}

    all_mers = list(mers.values())
    for from_mer in all_mers:
        for to_mer_name, bond_count in from_mer.bond_count.items():
            if bond_count > 0:
                if to_mer_name not in mers:
                    logger.warning(f"Mer {to_mer_name} not found in mers dictionary.")
                    continue
                to_mer_obj = mers[to_mer_name]
                affinity = bond_count / math.sqrt(len(from_mer.atoms) * len(to_mer_obj.atoms))
                weight = 1.0 / affinity
                adjacency_map[from_mer.name][to_mer_obj.name] = weight
                adjacency_map[to_mer_obj.name][from_mer.name] = weight

    # Remove Mers with no connections
    adjacency_map = {m: adj for m, adj in adjacency_map.items() if adj}
    return adjacency_map

def dijkstra(adjacency_map, source_mer):
    distances = {m: math.inf for m in adjacency_map}
    distances[source_mer] = 0.0
    predecessors = {m: None for m in adjacency_map}
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
                    predecessors[neighbor] = current

    return distances, predecessors

def generate_enhanced_pdb(weights, source_mer, original_file, mers):
    from io import StringIO
    output_buffer = StringIO()

    # Copy original header lines until first ATOM/HETATM
    with open(original_file, 'r') as reader:
        for line in reader:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                break
            output_buffer.write(line)

    output_buffer.write(f"REMARK Weights from source mer {source_mer} are stored in the temperature factor column.\n")
    output_buffer.write("REMARK Other fields such as atom ID, charge, and occupancy may be placeholders.\n\n")

    # Write out a single CA line per Mer with the distance in the B-factor
    for mer, weight in weights.items():
        if mer in mers and mers[mer].center_of_mass is not None:
            loc = mers[mer].center_of_mass
        else:
            loc = Location(0.0, 0.0, 0.0)

        try:
            residue_number = int(mer.split('-')[1].split('(')[0])
        except:
            residue_number = 0

        try:
            chain_id = mer.split('(')[1][0]
        except:
            chain_id = ' '

        output_buffer.write(
            "ATOM  {:>5d} {:<4s} {:>3s} {:1s}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:2s}\n".format(
                0,
                "CA",
                mer.split('-')[0],
                chain_id,
                residue_number,
                loc.x, loc.y, loc.z,
                1.0,
                weight,
                ""
            )
        )

    return output_buffer.getvalue()

def process_pdb_file(file_path):
    raw_lines = read_input_file(file_path)
    mers = parse_atoms(raw_lines)
    interactions = calculate_interactions(mers)
    adjacency_map = build_adjacency_map(mers)

    # Compute distances from every Mer to every other Mer
    all_distance_sums = {}
    best_source_mer = None
    min_total_distance = math.inf

    for mer in adjacency_map.keys():
        distances_source, _ = dijkstra(adjacency_map, mer)
        all_distance_sums[mer] = distances_source
        # Sum up reachable distances to see if this is the best source
        total_distance = sum(d for d in distances_source.values() if not math.isinf(d))
        if total_distance < min_total_distance:
            min_total_distance = total_distance
            best_source_mer = mer

    # For the best source Mer, we also create a simpler dictionary
    # that shows total_weight_sums from that best Mer
    distances_best = all_distance_sums[best_source_mer]
    total_weight_sums = {}
    for mer in distances_best:
        total_weight_sums[mer] = distances_best[mer]

    # Create an enhanced PDB for every possible Mer as source
    all_enhanced_pdbs = {}
    for mer, dist_map in all_distance_sums.items():
        enhanced_pdb = generate_enhanced_pdb(dist_map, mer, file_path, mers)
        all_enhanced_pdbs[mer] = enhanced_pdb

    return best_source_mer, mers, total_weight_sums, all_enhanced_pdbs, interactions, all_distance_sums
