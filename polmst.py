import math
import logging
import heapq
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
    """
    Return { mer_name: { neighbor_name: weight, ... }, ... }
    for all Mers, including possibly disconnected ones.
    """
    adjacency_map = {}
    for mer in mers.values():
        adjacency_map[mer.name] = {}

    for from_mer in mers.values():
        for to_mer_name, bond_count in from_mer.bond_count.items():
            if bond_count > 0 and to_mer_name in mers:
                to_mer = mers[to_mer_name]
                affinity = bond_count / math.sqrt(len(from_mer.atoms)*len(to_mer.atoms))
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
    """
    Returns a list of components, where each component is a set of Mer names.
    """
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
    """
    Build a sub‐map that has only the nodes in keep_set
    and edges among them.
    """
    new_map = {}
    for m in keep_set:
        if m in adjacency_map:
            # keep only neighbors also in keep_set
            sub_nbrs = {}
            for nbr, w in adjacency_map[m].items():
                if nbr in keep_set:
                    sub_nbrs[nbr] = w
            new_map[m] = sub_nbrs
    return new_map

def prune_interactions(interactions, sub_map):
    """
    Keep only those interactions whose from/to are in sub_map.
    """
    filtered = []
    for inter in interactions:
        if inter.from_mer in sub_map and inter.to_mer in sub_map:
            filtered.append(inter)
    return filtered

def prune_mers(mers, sub_map):
    """
    Keep only Mer objects in sub_map.
    """
    filtered = {}
    for m in sub_map:
        if m in mers:
            filtered[m] = mers[m]
    return filtered

def generate_enhanced_pdb(weights, source_mer, original_file, mers):
    from io import StringIO
    output_buffer = StringIO()

    with open(original_file, 'r') as reader:
        for line in reader:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                break
            output_buffer.write(line)

    output_buffer.write(f"REMARK Weights from source mer {source_mer} in TempFactor.\n\n")

    for mer, dist_val in weights.items():
        if mer in mers and mers[mer].center_of_mass is not None:
            loc = mers[mer].center_of_mass
        else:
            loc = Location(0,0,0)

        try:
            residue_number = int(mer.split('-')[1].split('(')[0])
        except:
            residue_number = 0
        try:
            chain_id = mer.split('(')[1][0]
        except:
            chain_id = ' '

        output_buffer.write(
            "ATOM  {:>5d} {:<4s} {:>3s} {:1s}{:>4d}    "
            "{:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:2s}\n".format(
                0, "CA", mer.split('-')[0],
                chain_id, residue_number,
                loc.x, loc.y, loc.z,
                1.0, dist_val, ""
            )
        )

    return output_buffer.getvalue()

def process_pdb_file(file_path):
    """
    1) Parse PDB -> build Mers
    2) Calculate interactions
    3) Build adjacency map
    4) Find all connected components -> pick the LARGEST
    5) Restrict adjacency_map, interactions, and mers to that largest component
    6) Compute best_source_mer in that largest component
    7) Do a final pass of Dijkstra for final distances
    """
    # 1) Parse
    raw_lines = read_input_file(file_path)
    mers = parse_atoms(raw_lines)
    # 2) Interactions
    interactions = calculate_interactions(mers)
    # 3) Adjacency
    adjacency_map = build_adjacency_map(mers)

    # 4) Find largest connected component
    comps = find_connected_components(adjacency_map)
    # pick the largest by node count
    largest_comp = max(comps, key=len)  # or handle tie breaks as you like

    # 5) Prune everything outside that largest component
    adjacency_map = prune_to_component(adjacency_map, largest_comp)
    interactions = prune_interactions(interactions, adjacency_map)
    mers = prune_mers(mers, adjacency_map)

    # 6) Find best_source_mer in that largest component
    #    by scanning all Mers in the largest comp
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

    # 7) The final distances from best_source_mer
    distances_best = all_distance_sums[best_source_mer]
    total_weight_sums = {m: distances_best[m] for m in distances_best}

    # Return everything
    return best_source_mer, mers, total_weight_sums, interactions, all_distance_sums
