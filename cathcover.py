import os
import re
import shutil
import argparse
import itertools
import subprocess
import numpy as np

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


def main(input_pdb, CATH_database, outfile, distance, coverage, tempfile):
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

    threshold_pairs = {}
    for line in lines[1:]:
        line_split = line.split()
        query_name = line_split[0]
        target_name = line_split[1]
        target_pdb = f"{results}/aligned_{query_name}_{target_name}.pdb"
        target_coords = CA_coords_from_pdb(target_pdb)
        target_nres = len(target_coords)
        qstart = int(line_split[2])
        tstart = int(line_split[4])
        cigar = line_split[8]
        pairs = cigar_to_pairs_from_seqs(query_nres, target_nres, cigar,
                                         qstart=qstart, tstart=tstart)
        threshold_pairs[target_name] = set()
        for pair in pairs:
            q_idx, t_idx = pair
            q_coord = query_coords[q_idx - 1]
            t_coord = target_coords[t_idx - 1]
            dist = np.linalg.norm(q_coord - t_coord)
            if dist <= distance:
                threshold_pairs[target_name].add(q_idx)

    # Calculate coverage
    cover_found = False
    for r in range(1, query_nres + 1):
        combs = itertools.combinations(threshold_pairs.keys(), r)
        for comb in combs:
            covered_residues = set()
            for target_name in comb:
                covered_residues.update(threshold_pairs[target_name])
            if len(covered_residues) / query_nres >= coverage:
                print(f"Coverage achieved with {r} CATH domains: {', '.join(comb)}")
                if outfile:
                    with open(outfile, 'a') as out_f:
                        out_f.write(f"{input_pdb}\t{','.join(comb)}\n")
                else:
                    print(f"{input_pdb}\t{','.join(comb)}")
                cover_found = True
                break
            if cover_found:
                break
        if cover_found:
            break

if __name__ == "__main__":
    args = parse_args()
    main(args.input_pdb, args.CATH_database, args.outfile,
         args.distance, args.coverage, args.tempfile)