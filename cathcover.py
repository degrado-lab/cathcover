import os
import re
import shutil
import argparse
import itertools
import subprocess
import numpy as np

from pathlib import Path
from typing import List, Sequence, Optional

def parse_args():
    parser = argparse.ArgumentParser(
        description="Determine the coverage of a protein by CATH domains."
    )
    parser.add_argument(
        'input_pdb', type=os.path.realpath, 
        help='Input PDB file path.'
    )
    parser.add_argument(
        'CATH_database', type=os.path.realpath, 
        help='Path to the CATH database file output by foldseek createdb.'
    )
    parser.add_argument(
        '--outfile', '-o', type=str, default=None,
        help='Output file to save the coverage results. '
        'If not specified, results are printed to stdout.'
    )
    parser.add_argument(
        '--pse-outfile', '-p', type=os.path.realpath, default=None,
        help='Path to output PSE file of aligned structures. If None, '
        'no PSE file will be created.'
    )
    parser.add_argument(
        '--distance', '-d', type=float, default=5.0,
        help='Distance threshold for coverage in Angstroms (default: 5.0).'
    )
    parser.add_argument(
        '--coverage', '-c', type=float, default=0.8,
        help='Proportion of residues that must be covered (default: 0.8).'
    )
    parser.add_argument(
        '--tempfile', '-t', type=os.path.realpath, default='/tmp',
        help='Temporary directory for Foldseek files (default: /tmp).'
    )
    return parser.parse_args()


def _as_resi_selector(resi_list: Sequence[int]) -> str:
    """Build a PyMOL selection string for a list of residue numbers."""
    if not resi_list:
        return "none"
    # Join as 'resi 1+2+3' (PyMOL handles this syntax)
    joined = "+".join(str(int(r)) for r in resi_list)
    return f"resi {joined}"


def _ensure_colors(cmd, n: int, names: Optional[List[str]] = None) -> List[str]:
    """
    Create n distinct color names in PyMOL and return their names.
    If `names` is given, it must have length n and will be used as-is.
    Otherwise, we create evenly spaced HSV colors and register them via set_color.
    """
    if names is not None:
        assert len(names) == n
        return names

    colors = []
    for i in range(n):
        # Evenly spaced hues; simple HSV->RGB conversion
        h = (i / max(1, n)) % 1.0
        s, v = 0.75, 0.95
        # HSV->RGB
        import colorsys
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        cname = f"group_color_{i+1}"
        cmd.set_color(cname, [r, g, b])
        colors.append(cname)
    return colors


def make_colored_pse(
    out_pse: str,
    main_pdb: str,
    other_pdbs: Sequence[str],
    main_resi_lists: Sequence[Sequence[int]],
    other_resi_lists: Sequence[Sequence[int]],
    object_names: Optional[Sequence[str]] = None,
    color_names: Optional[Sequence[str]] = None,
    show_as: str = "cartoon",
):
    """Create a colored PyMOL session (.pse) from the provided structures. 
    
    Parameters
    ----------
    out_pse : str
        Path to write the .pse session.
    main_pdb : str
        Path to primary PDB.
    other_pdbs : list[str]
        Paths to aligned PDBs (length N).
    main_resi_lists : list[list[int]]
        For each i, residues in main_pdb to color with the i-th color.
    other_resi_lists : list[list[int]]
        For each i, residues in other_pdbs[i] to color with the i-th color.
    object_names : optional list[str]
        Names to assign to loaded objects. If provided, length must be N+1 (main + each other).
        Defaults: ["main", "other_1", "other_2", ...]
    color_names : optional list[str]
        Custom PyMOL color names to use (or create). Must have length N if provided.
    show_as : str
        PyMOL representation for all objects (e.g., "cartoon", "sticks", "surface").
    """
    N = len(other_pdbs)
    if not (len(main_resi_lists) == len(other_resi_lists) == N):
        raise ValueError("Length mismatch: other_pdbs, main_resi_lists, other_resi_lists must be the same length")

    # Try to get a cmd handle whether we're in PyMOL or using pymol2
    try:
        # Case 1: Running *inside* PyMOL
        from pymol import cmd
        in_pymol = True
        pm = None
    except Exception:
        in_pymol = False
        pm = None

    if not in_pymol:
        # Case 2: Headless / normal Python using pymol2
        try:
            import pymol2  # type: ignore
        except Exception as e:
            raise RuntimeError(
                "Could not import PyMOL. Install PyMOL or run this script inside PyMOL.\n"
                "For headless usage: pip install pymol2"
            ) from e
        pm = pymol2.PyMOL()
        pm.start()
        from pymol import cmd  # now available

    try:
        cmd.reinitialize()  # clean session

        # Object naming
        if object_names is None:
            names = ["main"] + [f"other_{i+1}" for i in range(N)]
        else:
            if len(object_names) != N + 1:
                raise ValueError("object_names must have length N+1 (main + each other)")
            names = list(object_names)

        # Load structures
        cmd.load(main_pdb, names[0])
        for i, pdb in enumerate(other_pdbs):
            cmd.load(pdb, names[i + 1])

        # Set representations
        cmd.hide("everything")
        for nm in names:
            cmd.show_as(show_as, nm)

        # Prepare colors
        colors = _ensure_colors(cmd, N, list(color_names) if color_names else None)

        # Color per pair
        for i in range(N):
            main_sel = f"{names[0]} and ({_as_resi_selector(main_resi_lists[i])})"
            other_sel = f"{names[i+1]} and ({_as_resi_selector(other_resi_lists[i])})"
            # Name selections for convenience (optional)
            msel_name = f"{names[0]}_group_{i+1}"
            osel_name = f"{names[i+1]}_group_{i+1}"
            cmd.select(msel_name, main_sel)
            cmd.select(osel_name, other_sel)
            # Apply same color
            cmd.color(colors[i], msel_name)
            cmd.color(colors[i], osel_name)

        # A little prettifying (optional)
        cmd.bg_color("white")
        cmd.orient("all")

        # Save session
        out_pse = str(Path(out_pse))
        cmd.save(out_pse)
        print(f"Saved PyMOL session to: {out_pse}")

    finally:
        # If we created a headless PyMOL, shut it down
        if pm is not None:
            pm.stop()


def CA_coords_from_pdb(pdb_file):
    """
    Extract Cα atom coordinates from a PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    np.ndarray
        Array of shape (N, 3) with Cα coordinates.
    """
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append((x, y, z))
    return np.array(coords)


def cigar_to_pairs_from_seqs(query_nres, target_nres, cigar, qstart=1, tstart=1, one_based=True):
    """
    Build aligned residue index pairs from ungapped sequences + CIGAR.

    Parameters
    ----------
    query_nres : int
        Length of the ungapped query sequence.
    target_nres : int
        Length of the ungapped target sequence.
    cigar : str
        CIGAR string with operations M/I/D, e.g. "6M6I15M5D...".
        M: advance both (aligned pairs), I: advance query only, 
        D: advance target only.
    qstart : int, default=1
        1-based start position of the alignment in the query sequence 
        (Foldseek's qstart).
    tstart : int, default=1
        1-based start position of the alignment in the target sequence 
        (Foldseek's tstart).
    one_based : bool, default=True
        Output indices as 1-based (Foldseek style). If False, 
        return 0-based indices.

    Returns
    -------
    list[tuple[int,int]]
        List of (q_idx, t_idx) for aligned residue *pairs* (“M” runs only).

    Raises
    ------
    ValueError
        If the CIGAR would step past the sequence lengths.
    """
    # Parse CIGAR into (length, op) tuples
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MID])', cigar)]
    # Current positions (1-based inside the full sequences)
    qpos = qstart
    tpos = tstart
    pairs = []

    for length, op in ops:
        if op == 'M':
            # Emit length pairs and advance both
            qend = qpos + length - 1
            tend = tpos + length - 1
            if qend > query_nres or tend > target_nres:
                raise ValueError("CIGAR 'M' run exceeds sequence length.")
            if one_based:
                pairs.extend(zip(range(qpos, qend + 1), range(tpos, tend + 1)))
            else:
                pairs.extend(zip(range(qpos - 1, qend), range(tpos - 1, tend)))
            qpos += length
            tpos += length

        elif op == 'I':
            # Insertion in query relative to target: advance query only
            qend = qpos + length - 1
            if qend > query_nres:
                raise ValueError("CIGAR 'I' run exceeds query length.")
            qpos += length

        elif op == 'D':
            # Deletion in query (i.e., insertion in target): advance target only
            tend = tpos + length - 1
            if tend > target_nres:
                raise ValueError("CIGAR 'D' run exceeds target length.")
            tpos += length

        else:
            raise ValueError(f"Unsupported CIGAR op: {op}")

    return pairs


def main(input_pdb, CATH_database, outfile, pse_outfile, 
         distance, coverage, tempfile):
    """Main function to run the CATH coverage analysis.
    
    Parameters
    ----------
    input_pdb : str
        Path to the input PDB file.
    CATH_database : str
        Path to the CATH database file.
    outfile : str or None
        Path to the output file for coverage results. If None, results are 
        printed to stdout.
    pse_outfile : str of None
        Path to output PSE file of aligned structures. If None, no PSE file 
        will be created.
    distance : float
        Distance threshold for coverage in Angstroms.
    coverage : float
        Proportion of residues that must be covered.
    tempfile : str
        Temporary directory for Foldseek files.
    """
    # Read Cα coordinates from the input PDB
    query_coords = CA_coords_from_pdb(input_pdb)
    query_nres = len(query_coords)

    # Create Foldseek database for query PDB
    queries = f"{tempfile}/queries"
    os.makedirs(queries, exist_ok=True)
    shutil.copy(input_pdb, queries)
    queryDB = f"{tempfile}/queryDB"
    subprocess.run(['foldseek', 'createdb', queries, queryDB], check=True)

    # Create results directory in tempfile and change to it
    results = f"{tempfile}/foldseek_results_0"
    counter = 1
    while os.path.exists(results):
        results = '_'.join(results.split('_')[:-1]) + f"_{counter}"
        counter += 1
    os.makedirs(results, exist_ok=True)
    os.chdir(results)

    # Run Foldseek search
    subprocess.run(['foldseek', 'search', queryDB, CATH_database, 'aln_db',
                    tempfile, '-a', '--alignment-type', '2'], check=True)

    # Check if aln_db was created
    while not os.path.exists('aln_db.index') or \
    not os.path.exists('aln_db.dbtype'):
        time.sleep(1)

    # Output m8 file of aligned results
    format_output = ('query,target,qstart,qend,tstart,tend,qaln,taln,'
                     'cigar,lddtfull,u,t,alntmscore,qtmscore,ttmscore')
    subprocess.run(['foldseek', 'convertalis', queryDB, CATH_database, 'aln_db', 
                    'result.m8', '--format-mode', '4', '--format-output', 
                    format_output], check=True)

    # Output aligned PDB files
    subprocess.run(['foldseek', 'convertalis', queryDB, CATH_database, 'aln_db',
                    'aligned_', '--format-mode', '5'], check=True)

    # Read the m8 file and process results
    with open('result.m8', 'r') as f:
        lines = f.readlines()

    threshold_q_idxs = {}
    threshold_t_idxs = {}
    target_pdbs = {}
    for line in lines[1:]:
        line_split = line.split()
        query_name = line_split[0]
        target_name = line_split[1]
        target_pdb = f"{results}/aligned_{query_name}_{target_name}.pdb"
        target_pdbs[target_name] = target_pdb
        target_coords = CA_coords_from_pdb(target_pdb)
        target_nres = len(target_coords)
        qstart = int(line_split[2])
        tstart = int(line_split[4])
        cigar = line_split[8]
        pairs = cigar_to_pairs_from_seqs(query_nres, target_nres, cigar,
                                         qstart=qstart, tstart=tstart)
        threshold_q_idxs[target_name] = []
        threshold_t_idxs[target_name] = []
        for pair in pairs:
            q_idx, t_idx = pair
            q_coord = query_coords[q_idx - 1]
            t_coord = target_coords[t_idx - 1]
            dist = np.linalg.norm(q_coord - t_coord)
            if dist <= distance:
                threshold_q_idxs[target_name].append(q_idx)
                threshold_t_idxs[target_name].append(t_idx)

    # Calculate coverage
    cover_found = False
    for r in range(1, query_nres + 1):
        combs = itertools.combinations(threshold_q_idxs.keys(), r)
        for comb in combs:
            prev_covered_residues_set = set()
            covered_residues_set = set()
            covered_residues = []
            matching_residues = []
            coverages = []
            for i, target_name in enumerate(comb):
                assert len(threshold_t_idxs[target_name]) == \
                    len(threshold_q_idxs[target_name])
                prev_covered_residues_set = covered_residues_set.copy()
                covered_residues_set.update(threshold_q_idxs[target_name])
                if i == 0:
                    covered_residues.append(list(covered_residues_set))
                else:
                    covered_residues.append(list(covered_residues_set.difference(
                        prev_covered_residues_set
                    )))
                matching_residues.append([
                    threshold_t_idxs[target_name][threshold_q_idxs[target_name].index(q_idx)] 
                    for q_idx in covered_residues[i]
                ])
                coverages.append(len(covered_residues[i]) / query_nres)
            percent_covered = sum(coverages)
            if percent_covered >= coverage:
                print((f"Coverage of {percent_covered:.2%} "
                       f"achieved with {r} CATH domains: {', '.join(comb)}"))
                print("Coverage details: " + ", ".join(f"{c:.2%}" 
                                                       for c in coverages))
                if outfile:
                    with open(outfile, 'a') as out_f:
                        out_f.write(f"{input_pdb}\t{','.join(comb)}\n")
                else:
                    print(f"{input_pdb}\t{','.join(comb)}")
                if pse_outfile:
                    match_pdbs = [target_pdbs[name] for name in comb]
                    # for cr, mr in zip(covered_residues, matching_residues):
                    #     print(cr, '\n', mr, '\n')
                    color_names = ['marine', 'firebrick', 'splitpea', 
                                   'gold', 'orange', 'deepviolet', 
                                   'aquamarine', 'sand', 'olive', 'salmon']
                    if len(match_pdbs) < len(color_names):
                        make_colored_pse(pse_outfile, input_pdb, match_pdbs, 
                                         covered_residues, matching_residues, 
                                         color_names=color_names[:len(match_pdbs)])
                    else:
                        make_colored_pse(pse_outfile, input_pdb, match_pdbs, 
                                         covered_residues, matching_residues)
                cover_found = True
                break
            if cover_found:
                break
        if cover_found:
            break

if __name__ == "__main__":
    args = parse_args()
    main(args.input_pdb, args.CATH_database, args.outfile,
         args.pse_outfile, args.distance, args.coverage, args.tempfile)