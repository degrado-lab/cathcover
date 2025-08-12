# cathcover
Determine the number of CATH domains required to greedily cover a PDB structure.

## Description
This repository follows the method used by [Ingraham et al.](https://www.nature.com/articles/s41586-023-06728-8) to determine the novelty of protein folds by determining the number of CATH domains required to greedily cover 80% of a fold at 5 Angstroms RMSD per cover. Optional command line arguments allow the user to adjust the 80% or 5-Angstrom thresholds.

## Installation
To use this repository, one must install NumPy into the local environment:

```bash
$ pip install numpy
```

Alternatively, one may create an environment with NumPy:

```bash
$ conda create --name env_numpy numpy
```

One must then install the Foldseek binaries:

```bash
$ cd
$ wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz
```

Then add the Foldseek binary folder to the `$PATH` environment variable by adding the following line to your `.bashrc`:

```bash
export PATH="$HOME/foldseek/bin:$PATH"
```

Remember to source the `.bashrc` afterwards to update your `$PATH`!

Lastly, a copy of the CATHdb S40 database should be downloaded as follows:

```bash
$ cd
$ wget ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/non-redundant-data-sets/cath-dataset-nonredundant-S40.pdb.tgz
$ tar -xvzf cath-dataset-nonredundant-S40.pdb.tgz
$ for f in $HOME/dompdb; do mv -- "$f" "$f.pdb"; done
```

## Usage
To use the software, simply provide the `cathcover.py` script with the path to the PDB file for which you wish to calculate the coverage and the path of the CATHdb folder you downloaded:

```bash
$ python cathcover/cathcover.py mypdb.pdb $HOME/dompdb
```

The script will run and return the number of CATH structures required to greedily cover the PDB file of interest at a 5-Angstrom RMSD threshold. By default, these are printed to the display, although the `-o`, or `--outfile`, argument may be used to specify an output file for this information. The `-r`, or `--rmsd`, argument may be used to adjust the per-cover RMSD threshold from a default of 5 Angstroms, and the `-c`, or `--coverage`, argument may be used to adjust the proportion of residues that must be covered from a default of 0.8.
