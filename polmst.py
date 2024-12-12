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
    # CHANGED: Incorporate logic from Java code.
    # We first identify bonds by checking COM distances and then individual atom pairs.
    # After that, we list interactions.
    mer_list = list(mers.values())

    for i in range(len(mer_list)):
        src = mer_list[i]
        for j in range(len(mer_list)):
            if i == j:
                continue
            dst = mer_list[j]

            # Check COM distance threshold of 10.0 Å
            if src.center_of_mass.distance_to(dst.center_of_mass) <= 10.0:
                # Check atoms pairwise
                bonded = False
                for src_atom in src.atoms:
                    for dst_atom in dst.atoms:
                        if src_atom.location.distance_to(dst_atom.location) < 4.5:
                            src.add_bond(dst)
                            dst.add_bond(src)
                            bonded = True
                            # Once we have found a bond, no need to check further atom pairs
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
    # CHANGED: Ensure symmetry as in the Java code.
    adjacency_map = {}
    for mer in mers.values():
        adjacency_map[mer] = {}

    all_mers = list(mers.values())
    for from_mer in all_mers:
        # Assume `from_mer.bond_count` is similar to the Java version: 
        # a dict mapping another mer to an integer bond count.
        for to_mer, bond_count in from_mer.bond_count.items():
            if bond_count > 0:
                affinity = bond_count / math.sqrt(len(from_mer.atoms) * len(to_mer.atoms))
                weight = 1.0 / affinity
                adjacency_map[from_mer][to_mer] = weight
                # Ensure symmetry
                if to_mer not in adjacency_map:
                    adjacency_map[to_mer] = {}
                adjacency_map[to_mer][from_mer] = weight

    # Remove any empty entries if needed
    adjacency_map = {m: adj for m, adj in adjacency_map.items() if adj}
    return adjacency_map

def dijkstra(adjacency_map, source_mer):
    distances = {m: math.inf for m in adjacency_map}
    distances[source_mer] = 0.0
    visited = set()
    unvisited = set(adjacency_map.keys())

    while unvisited:
        # Pick the mer with the smallest known distance that is still unvisited
        current = min(unvisited, key=lambda m: distances[m])
        unvisited.remove(current)
        visited.add(current)

        # If the smallest distance is infinity, remaining nodes are unreachable
        if distances[current] == math.inf:
            break

        # Update distances to neighbors
        for neighbor, weight in adjacency_map[current].items():
            if neighbor not in visited:
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

import os
import math

import os
import math
from io import StringIO

def generate_enhanced_pdb(distances, source_mer, original_file):
    """
    Generate a PDB-formatted string annotated with distances from a given source mer.
    Distances are stored in the B-factor (temperature factor) column of the ATOM records.

    Parameters:
        distances (dict): A dictionary mapping Mer objects to their shortest distance from source_mer.
        source_mer (Mer): The reference Mer object from which distances were calculated.
        original_file (str or Path): The path to the original PDB file.

    Returns:
        str: The resulting PDB-formatted string with distances annotated.
    """

    # Use a string buffer to build the output
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
    output_buffer.write("REMARK Distances from source mer {} are stored in the temperature factor column.\n".format(source_mer.name))
    output_buffer.write("REMARK Other fields such as atom ID, charge, and occupancy may be placeholders.\n\n")

    # Write out the distances as ATOM lines
    # Each Mer is represented by a single CA line with the distance as the B-factor.
    for mer, distance in distances.items():
        # Retrieve the center of mass location for coordinates (fallback to 0,0,0 if None)
        loc = mer.center_of_mass if mer.center_of_mass is not None else type('Location', (), {'x':0.0, 'y':0.0, 'z':0.0})()

        output_buffer.write(
            "ATOM  {:>5d} {:<4s} {:>3s} {:1s}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:2s}\n".format(
                0,                          # Atom serial number placeholder
                "CA",                       # Atom name
                mer.type,                   # Residue name
                mer.chain,                  # Chain ID
                mer.mer_id,                 # Residue sequence number
                loc.x, loc.y, loc.z,        # Coordinates (from center of mass)
                1.0,                        # Occupancy
                distance,                   # B-factor = distance
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
    best_source_mer = find_best_source_mer(mers.values(), adjacency_map)

    # Generate enhanced PDB for every mer as the source
    all_enhanced_pdbs = {}
    all_distances = {}
    for source_mer in mers.values():
        distances = dijkstra(adjacency_map, source_mer)
        enhanced_pdb = generate_enhanced_pdb(distances, source_mer, file_path)
        all_enhanced_pdbs[source_mer.name] = enhanced_pdb
        all_distances[source_mer.name] = distances

    shortest_paths = all_distances[best_source_mer.name]
    return best_source_mer.name, mers, shortest_paths, all_enhanced_pdbs
