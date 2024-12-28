# /Users/m3_air/24-25 fall/Bitirme/bito-pro-main/polmst.py

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

    # Now build the list of interactions from bonded Mers
    interactions = []
    for i in range(len(mer_list)):
        for j in range(i+1, len(mer_list)):
            src = mer_list[i]
            dst = mer_list[j]
            if src.is_bonded_to(dst):
                inter = Interaction(src, dst)
                interactions.append(inter)

    return interactions

def build_adjacency_map(mers):
    """
    Builds an adjacency map where each Mer is a node, and edges represent interactions
    weighted by their calculated weights. Ensures that the adjacency map is symmetric.
    """
    adjacency_map = {}
    for mer in mers.values():
        adjacency_map[mer.name] = {}  # Initialize adjacency list for each Mer

    all_mers = list(mers.values())
    for from_mer in all_mers:
        for to_mer_name, bond_count in from_mer.bond_count.items():
            if bond_count > 0:
                if to_mer_name not in mers:
                    # Log warning and skip if to_mer not found
                    logger.warning(f"Mer {to_mer_name} not found in mers dictionary.")
                    continue
                to_mer_obj = mers[to_mer_name]  # Look up the Mer object using its name
                # Calculate affinity
                affinity = bond_count / math.sqrt(len(from_mer.atoms) * len(to_mer_obj.atoms))
                # Calculate weight as inverse of affinity
                weight = 1.0 / affinity
                adjacency_map[from_mer.name][to_mer_obj.name] = weight
                # Ensure symmetry in the adjacency map
                adjacency_map[to_mer_obj.name][from_mer.name] = weight

    # Remove any Mers with no connections
    adjacency_map = {m: adj for m, adj in adjacency_map.items() if adj}
    return adjacency_map


def dijkstra(adjacency_map, source_mer):
    """
    Implements Dijkstra's algorithm to compute the shortest paths from the source Mer
    to all other Mers in the adjacency map.

    Parameters:
        adjacency_map (dict): The adjacency map with weights.
        source_mer (str): The name of the source Mer.

    Returns:
        tuple: A dictionary of shortest distances and a dictionary of predecessors.
    """
    distances = {m: math.inf for m in adjacency_map}
    distances[source_mer] = 0.0
    predecessors = {m: None for m in adjacency_map}
    visited = set()
    unvisited = set(adjacency_map.keys())

    while unvisited:
        # Select the unvisited Mer with the smallest known distance
        current = min(unvisited, key=lambda m: distances[m])
        unvisited.remove(current)
        visited.add(current)

        # If the smallest distance is infinity, remaining Mers are unreachable
        if distances[current] == math.inf:
            break

        # Update distances to neighboring Mers
        for neighbor, weight in adjacency_map[current].items():
            if neighbor not in visited:
                new_dist = distances[current] + weight
                if new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    predecessors[neighbor] = current

    return distances, predecessors

def find_best_source_mer(mers, adjacency_map):
    """
    Identifies the best source Mer based on the smallest total weight to all other Mers.

    Parameters:
        mers (dict): Dictionary of Mer objects.
        adjacency_map (dict): The adjacency map with weights.

    Returns:
        str: The name of the best source Mer.
    """
    best_mer = None
    min_total_distance = math.inf
    for source_mer in adjacency_map.keys():
        dist, _ = dijkstra(adjacency_map, source_mer)
        total_distance = sum(d for d in dist.values() if not math.isinf(d))
        if total_distance < min_total_distance:
            min_total_distance = total_distance
            best_mer = source_mer
    return best_mer

def reconstruct_path(predecessors, target_mer):
    """
    Reconstructs the shortest path from the source Mer to the target Mer.

    Parameters:
        predecessors (dict): Dictionary mapping each Mer to its predecessor in the path.
        target_mer (str): The name of the target Mer.

    Returns:
        list: The list of Mers representing the shortest path.
    """
    path = []
    current = target_mer
    while current is not None:
        path.append(current)
        current = predecessors[current]
    path.reverse()
    return path

def sum_weights_along_path(mers, interactions, path):
    """
    Sums the weights of interactions along a given path of Mers.

    Parameters:
        mers (dict): Dictionary of Mer objects.
        interactions (list): List of Interaction objects.
        path (list): The list of Mers representing the path.

    Returns:
        float: The total weight sum along the path.
    """
    total_weight = 0.0
    for i in range(len(path) -1 ):
        from_mer = path[i]
        to_mer = path[i+1]
        # Find the interaction between from_mer and to_mer
        for interaction in interactions:
            if (interaction.from_mer == from_mer and interaction.to_mer == to_mer) or \
               (interaction.from_mer == to_mer and interaction.to_mer == from_mer):
                total_weight += interaction.weight
                break
    return total_weight

def generate_enhanced_pdb(weights, source_mer, original_file, mers):
    """
    Generates a PDB-formatted string annotated with weights from a given source Mer.
    Weights are stored in the B-factor (temperature factor) column of the ATOM records.

    Parameters:
        weights (dict): A dictionary mapping Mer names to their shortest weight from source_mer.
        source_mer (str): The reference Mer name from which weights were calculated.
        original_file (str or Path): The path to the original PDB file.
        mers (dict): Dictionary of Mer objects with mer_name as key.

    Returns:
        str: The resulting PDB-formatted string with weights annotated.
    """

    # Use a string buffer to build the output
    from io import StringIO
    output_buffer = StringIO()

    # Read the original file
    with open(original_file, 'r') as reader:
        # Copy headers from the original file until the first ATOM/HETATM line
        for line in reader:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Once we reach an ATOM/HETATM line, stop copying headers
                # We don't need to process the original ATOM/HETATM lines for this output,
                # so we break here.
                break
            output_buffer.write(line)

    # Write explanatory headers
    output_buffer.write("REMARK Weights from source mer {} are stored in the temperature factor column.\n".format(source_mer))
    output_buffer.write("REMARK Other fields such as atom ID, charge, and occupancy may be placeholders.\n\n")

    # Write out the weights as ATOM lines
    # Each Mer is represented by a single CA line with the weight as the B-factor.
    for mer, weight in weights.items():
        # Retrieve the center of mass location for coordinates (fallback to 0,0,0 if None)
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
                0,                          # Atom serial number placeholder
                "CA",                       # Atom name
                mer.split('-')[0],          # Residue name
                chain_id,       # Chain ID
                residue_number,  # Residue sequence number
                loc.x, loc.y, loc.z,        # Coordinates (from center of mass)
                1.0,                        # Occupancy
                weight,                     # B-factor = weight
                ""                          # Element symbol
            )
        )

    # Return the entire PDB content as a string
    return output_buffer.getvalue()



def process_pdb_file(file_path):
    raw_lines = read_input_file(file_path)
    mers = parse_atoms(raw_lines)
    interactions = calculate_interactions(mers)
    adjacency_map = build_adjacency_map(mers)

    # Find the best source mer
    best_source_mer = find_best_source_mer(mers, adjacency_map)
    if best_source_mer is None:
        raise ValueError("No valid source mer found.")

    # Perform Dijkstra from the best source mer
    distances, predecessors = dijkstra(adjacency_map, best_source_mer)

    # For each mer, reconstruct the path and sum weights
    total_weight_sums = {}
    for mer in adjacency_map.keys():
        if mer == best_source_mer:
            total_weight_sums[mer] = 0.0  # Source mer has sum 0
            continue
        if math.isinf(distances[mer]):
            total_weight_sums[mer] = math.inf  # Unreachable
            continue
        path = reconstruct_path(predecessors, mer)
        total_weight = sum_weights_along_path(mers, interactions, path)
        total_weight_sums[mer] = total_weight

    # Generate enhanced PDBs
    all_enhanced_pdbs = {}
    for source_mer in adjacency_map.keys():
        distances_source, _ = dijkstra(adjacency_map, source_mer)
        enhanced_pdb = generate_enhanced_pdb(distances_source, source_mer, file_path, mers)
        all_enhanced_pdbs[source_mer] = enhanced_pdb

    # Return the best source mer's name, mers, total_weight_sums, all_enhanced_pdbs, interactions
    return best_source_mer, mers, total_weight_sums, all_enhanced_pdbs, interactions
