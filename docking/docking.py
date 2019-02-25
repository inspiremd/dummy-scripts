#!/usr/bin/env python3
import os
from argparse import ArgumentParser
import tempfile
from openeye import oechem, oedocking

def dock_smiles(smiles_filename, pdb_filename, output_dir, **kwargs):
    """
    Dummy function to stand in for Docking module of INSPIRE MD/ML combined
    workflow. The module will eventually both dock input ligands into provided
    protein (PDB) structure and produce system prepared for MD smiulation 
    using tLeap (AMBER).

    Parameters
    ----------
    smiles_filename: str_like
        Path to the input file contaning one or more SMILES strings.
    pdb_filename: str_like
        Path to PDB file into which to dock componds represented by the
        input SMILES.
    output_dir: str_like
        Path where directories will be created to store output of docking and 
        parameterization.
        
    Returns
    -------
    array
        List of paths within output_dir populated with 
    """

    tempfile.tempdir = output_dir

    created_systems = []

    with open(smiles_filename) as input_smiles:
        
        smiles = input_smiles.readline().strip()

        ligand_dir = tempfile.mkdtemp()

        with open(os.path.join(ligand_dir, 'smiles.txt', 'w') as output_smiles:
            print(smiles, file=output_smiles)

        with open(os.path.join(ligand_dir, 'complex.prmtop', 'w') as output_top:
            print("This would be an AMBER topology file", file=output_top)
        
        with open(os.path.join(ligand_dir, 'complex.inpcrd', 'w') as output_crd:
            print("This would be an AMBER input coordinates files", file=output_crd)
        
        with open(os.path.join(ligand_dir, 'complex.pdb', 'w') as output_pdb:
            print("This would be an PDB input coordinates files", file=output_pdb)

        with open(os.path.join(ligand_dir, 'stripped_complex.prmtop', 'w') as output_top:
            print("This would be an AMBER topology file for system without solvent", 
                  file=output_top)
        
        with open(os.path.join(ligand_dir, 'stripped_complex.pdb', 'w') as output_pdb:
            print("This would be an PDB input coordinates files for system without solvent", 
                  file=output_pdb)

        ligand_specific_dir = os.path.split(ligand_dir)[-1]
        created_systems.append(ligand_specific_dir)


    return created_systems

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-s", "--smiles", dest="smiles_filename",
                        help="File containing lists of SMILES strings")
    parser.add_argument("-p", "--pdb", dest="pdb_filename",
                        help="File containing PDB to dock SMILES into")
    parser.add_argument("-o", "--output_dir", dest="output_dir", 
                        help="Directory into which the simulation initiation "
                        "files will be output")

    args = parser.parse_args()

    dock_smiles(args.smiles_filename, args.pdb_filename, args.output_dir)
