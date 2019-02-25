#!/usr/bin/env python3
import os
from argparse import ArgumentParser
import numpy as np
from scipy.stats import halfnorm
# TODO: Add relevant OpenMM imports

def dg_from_dist(n_samples=25, offset=-10.0):
    """
    Generate some fake results for MMGBSA then return mean and
    standard error as an uncertainty estimate.

    Parameters
    ----------
    n_samples: integer
        Number of samples to draw
    offset: float
        Offset from 0 - generally MMGBSA results are more negative than 
        experimental ones
    """

    target = -1 * halfnorm.rvs(size=1)[0] + offset

    sample = np.random.randn(n_samples) + target

    return np.mean(sample), np.sd(sample)/np.sqrt(n_samples)

def run_mmgbsa(input_dir, **kwargs):
    """
    Dummy function to stand in MD run and MMGBSA analysis module of INSPIRE MD/ML 
    combined workflow. The module will eventually both dock input ligands into 
    provided protein (PDB) structure and produce system prepared for MD 
    smiulation using tLeap (AMBER).

    Parameters
    ----------
    input_dir: str_like
        Path containing input coordinates and topology.
    """

    # NOTE: In final version we may well make multiple directories for 
    #       replica runs.

    run_dir = os.path.join(input_dir, 'run1')
    os.makdirs(run_dir)
    
    analysis_dir = os.path.join(input_dir, 'stats')
    os.makdirs(analysis_dir)

    with open(smiles_filename) as input_smiles:
        
        smiles = input_smiles.readline().strip()

        ligand_dir = tempfile.mkdtemp()

        with open(os.path.join(run_dir, 'output.dcd', 'w') as output_traj:
            print("Trajectory", file=output_traj)

        with open(os.path.join(run_dir, 'stripped_output.dcd', 'w') as output_traj:
            print("Solvent stripped trajectory", file=output_traj)

        with open(os.path.join(run_dir, 'mmgbsa.dat', 'w') as output_mmgbsa:
            print("This would be an AMBER topology file", file=output_mmgbsa)

        # TODO: Generate values from a distribution

        dg, err = dg_from_dist()
        
        with open(os.path.join(analysis_dir, 'summary_mmgbsa.dat', 'w') as output_stats:
            print("{:.2f}\t{:2f}".format(dg, err), file=output_stats)

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", dest="input_dir",
                        help="Directory contining simulation inputs (topology and coordinates)")

    args = parser.parse_args()

    run_mmgbsa(args.input_dir)
