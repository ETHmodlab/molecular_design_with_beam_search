# Copyright (c) 2021 ETH Zurich

import sys, os
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def get_mol(smi):
    return Chem.MolFromSmiles(smi)

def get_Morgan(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
    return fp

def extract_murcko_scaffolds(mols, verbose=True):
    """ Extract Bemis-Murcko scaffolds from a smile string.

    :param mols: molecule data set in rdkit mol format.
    :return: smiles string of a scaffold and a framework.
    """
    scaf = []
    scaf_unique = []
    generic_scaf = []
    generic_scaf_unique = []
    start = time.time()
    for mol in mols:
        if mol is None:
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            fw = MurckoScaffold.MakeScaffoldGeneric(core)
            scaf.append(Chem.MolToSmiles(core, isomericSmiles=True))
            generic_scaf.append(Chem.MolToSmiles(fw, isomericSmiles=True))
        except ValueError as e:
            print(e)
            scaf.append(['error'])
            generic_scaf.append(['error'])
    if verbose:
        print('Extracted', len(scaf), 'scaffolds in', time.time() - start, 'seconds.')
    return scaf, generic_scaf
        
      