# cathcover
Determine the number of CATH domains required to greedily cover a PDB structure.

## Description
This repository follows the method used by [Ingraham et al.](https://www.nature.com/articles/s41586-023-06728-8) to determine the novelty of a query protein fold by determining the number of CATH domains required to greedily cover 80% of a fold with residues that lie within a CA-CA distance of 5 Angstroms of the corresponding residue in the query. Optional command line arguments allow the user to adjust the 80% or 5-Angstrom thresholds.

## Installation
To use this repository, one must install NumPy into the local environment:

```bash
$ pip install numpy
```

Alternatively, one may create an environment with NumPy:

```bash
$ conda create --name env_cathcover numpy
```

Optionally, one can also ensure that PyMOL is installed into one's environment to enable output of PSEs with the aligned CATH domains (vide infra):

```bash
pip install numpy pymol-open-source
```

Or, via conda:

```bash
conda create --name env_cathcover numpy pymol-open-source
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

To accelerate future Foldseek searches, it is helpful to create a Foldseek database object:

```bash
$ mkdir CATHdb
$ cd CATHdb
$ foldseek createdb $HOME/dompdb CATHdb
```

## Usage
To use the software, simply provide the `cathcover.py` script with the path to the PDB file for which you wish to calculate the coverage and the path of the CATHdb database file you created:

```bash
$ python cathcover/cathcover.py mypdb.pdb $HOME/CATHdb/CATHdb
```

The script will run and return the number of CATH structures required to greedily cover the PDB file of interest at a 5-Angstrom CA-CA distance threshold. By default, these are printed to the display, although the `-o`, or `--outfile`, argument may be used to specify an output file for this information. The `-d`, or `--distance`, argument may be used to adjust the per-cover distance threshold from a default of 5 Angstroms, and the `-c`, or `--coverage`, argument may be used to adjust the proportion of residues that must be covered from a default of 0.8. The `-t`, or `--tempfile`, argument specifies the temporary directory in which Foldseek files will be output, with `/tmp` as the default. Lastly, the `-p`, or `--pse-outfile`, argument allows one to output a PSE file with the aligned structures of the covering CATH domains shown alongside the query protein and the matching residues colored in the same color. This option requires PyMOL to be in the environment.
