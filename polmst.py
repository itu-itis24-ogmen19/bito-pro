import math
import logging
from io import StringIO

from models import Atom, Mer, Location, Interaction

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# ────────────────────────────────────────────────────────────────
#  PDB parsing helpers
# ────────────────────────────────────────────────────────────────
def read_input_file(file_path):
    with open(file_path, "r") as f:
        return f.readlines()


def parse_atoms(raw_lines):
    """Create Mer objects and attach Atom objects to them."""
    mers = {}
    for line in raw_lines:
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_atom_line(line)
            if atom:
                mers.setdefault(atom.mer.name, atom.mer).add_atom(atom)

    for mer in mers.values():
        mer.calc_com()
    return mers


def parse_atom_line(line):
    """Return an Atom or None (if line malformed)."""
    try:
        if len(line) < 66:
            return None
        atom_id = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resnum = int(line[22:26])
        x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
        temp = float(line[60:66])

        loc = Location(x, y, z)
        mer_name = f"{resname}-{resnum}({chain})"
        mer = Mer(mer_name, resname, resnum, chain)
        return Atom(atom_id, name, resname, "", loc, temp, mer)
    except Exception:
        return None


# ────────────────────────────────────────────────────────────────
#  Graph helpers
# ────────────────────────────────────────────────────────────────
def calculate_interactions(mers):
    """Mark bonds (<4.5 Å) and create Interaction objects."""
    mer_list = list(mers.values())

    # register bonds
    for i, src in enumerate(mer_list):
        for j, dst in enumerate(mer_list):
            if i == j:
                continue
            bonded = any(
                s_atom.location.distance_to(d_atom.location) < 4.5
                for s_atom in src.atoms
                for d_atom in dst.atoms
            )
            if bonded:
                src.add_bond(dst)
                dst.add_bond(src)

    # build Interaction list
    interactions = []
    for i in range(len(mer_list)):
        for j in range(i + 1, len(mer_list)):
            a, b = mer_list[i], mer_list[j]
            if a.is_bonded_to(b):
                interactions.append(Interaction(a, b))
    return interactions


def build_adjacency_map(mers):
    """Return {mer: {nbr: weight}} for every Mer."""
    adj = {m.name: {} for m in mers.values()}
    for a in mers.values():
        for b_name, bond_count in a.bond_count.items():
            if bond_count and b_name in mers:
                b = mers[b_name]
                affinity = bond_count / math.sqrt(len(a.atoms) * len(b.atoms))
                weight = 1.0 / affinity
                adj[a.name][b_name] = weight
                adj[b_name][a.name] = weight
    return adj


def dijkstra(adj, src):
    """Return distance map from src."""
    dist = {m: math.inf for m in adj}
    dist[src] = 0.0
    visited = set()

    while (unvisited := {m for m in adj if m not in visited}):
        current = min(unvisited, key=dist.get)
        visited.add(current)
        if dist[current] == math.inf:
            break
        for nbr, w in adj[current].items():
            if nbr not in visited and dist[current] + w < dist[nbr]:
                dist[nbr] = dist[current] + w
    return dist


def find_connected_components(adj):
    comps, seen = [], set()
    for node in adj:
        if node in seen:
            continue
        stack, comp = [node], set()
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.add(n)
            stack.extend(adj[n])
        comps.append(comp)
    return comps


def prune_to_component(adj, keep):
    return {n: {k: w for k, w in nbrs.items() if k in keep} for n, nbrs in adj.items() if n in keep}


def prune_interactions(interactions, keep_map):
    return [i for i in interactions if i.from_mer in keep_map and i.to_mer in keep_map]


def prune_mers(mers, keep_map):
    return {m: mers[m] for m in keep_map if m in mers}


# ────────────────────────────────────────────────────────────────
#  PDB writer (keeps all atoms, updates B-factor)
# ────────────────────────────────────────────────────────────────
def generate_enhanced_pdb(weights, source_mer, original_file, _mers):
    """
    Copy every original atom line.  Put distance/weight for the atom’s Mer
    into the B-factor (TempFactor) column.  Add one REMARK before coords.
    """
    out = StringIO()
    with open(original_file, "r") as pdb:
        remark_done = False
        for line in pdb:
            is_atom = line.startswith(("ATOM", "HETATM"))
            if is_atom and not remark_done:
                out.write(
                    f"REMARK 999  Distances from {source_mer} are stored in the B-factor column\n"
                )
                remark_done = True

            if is_atom and len(line) >= 66:
                res = line[17:20].strip()
                chain = line[21].strip()
                num = line[22:26].strip()
                mer_name = f"{res}-{num}({chain})"
                bf = weights.get(mer_name)
                if bf is not None and math.isfinite(bf):
                    bf = min(bf, 999.99)
                    line = f"{line[:60]}{bf:6.2f}{line[66:]}"
            out.write(line)
    return out.getvalue()


# ────────────────────────────────────────────────────────────────
#  High-level processing pipeline
# ────────────────────────────────────────────────────────────────
def process_pdb_file(file_path):
    """
    Returns:
        best_source_mer, mers, total_weight_sums, interactions,
        all_distance_sums, islands_removed (bool)
    """
    raw = read_input_file(file_path)
    mers = parse_atoms(raw)

    interactions = calculate_interactions(mers)
    adj = build_adjacency_map(mers)

    comps = find_connected_components(adj)
    largest = max(comps, key=len)
    islands_removed = len(comps) > 1

    adj = prune_to_component(adj, largest)
    interactions = prune_interactions(interactions, adj)
    mers = prune_mers(mers, adj)

    all_dist, best_src, min_sum = {}, None, math.inf
    for mer in largest:
        d = dijkstra(adj, mer)
        all_dist[mer] = d
        s = sum(v for v in d.values() if math.isfinite(v))
        if s < min_sum:
            best_src, min_sum = mer, s

    total_weight_sums = dict(all_dist[best_src])
    return best_src, mers, total_weight_sums, interactions, all_dist, islands_removed
