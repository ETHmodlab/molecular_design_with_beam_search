# Copyright (c) 2021 ETH Zurich

import os, sys
import argparse
import time
import configparser
import warnings

from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
from rdkit.Chem import Draw

from chembl_webresource_client.new_client import new_client

sys.path.append('../src/')
from python import helper as hp

parser = argparse.ArgumentParser(description='Get the final ranking')
parser.add_argument('-c','--configfile', type=str, help='Path to config file', required=True)
parser.add_argument('-v','--verbose', type=bool, help='Verbose', required=True)


def is_novel_chembl(smi):
    """
    Return True if the SMILES smi 
    was found in ChEMBL, False if not
    """
    similarity = new_client.similarity
    res = similarity.filter(smiles=smi, similarity=100)
    if res:
        return False
    else:
        return True
    
       
if __name__ == '__main__':
    
    start = time.time()
    
    ####################################
    # get back parameters
    args = vars(parser.parse_args())
    
    verbose = args['verbose']
    configfile = args['configfile']
    config = configparser.ConfigParser()
    config.read(configfile)
    exp_name = configfile.split('/')[-1].replace('.ini','')
    
    if verbose: print('\nSTART FINAL RANKING')
    ####################################
    
    
    
    
    ####################################
    # Path to save the ranking
    save_path = f'results/{exp_name}/final_ranking/'
    os.makedirs(save_path, exist_ok=True)
    ####################################

    
    
    
    ####################################
    # start iterating over the files
    dir_experiment = f'results/{exp_name}/sampling/'
    # first fine-tuning epoch from which
    # we start considering molecules for the
    # final rank
    from_epoch = int(config['BEAM']['from_epoch'])
    
    d_smi_score = {}
    for filename in os.listdir(dir_experiment):
        if filename.endswith('.pkl'):
            epoch = int(filename.split('_')[0])
            # we only consider epoch from the
            # epoch defined as from_epoch
            if epoch >= from_epoch:
                _d = hp.load_obj(f'{dir_experiment}{filename}')
                # we update the general dict.
                # if we have a smi appearing twice,
                # we keep the best score
                for k,v in _d.items():
                    if k in d_smi_score:
                        if v>d_smi_score[k]:
                            d_smi_score[k] = v
                    else:
                        d_smi_score[k] = v
                        
    # sort the results
    sorted_score = sorted(d_smi_score.items(), key=lambda x: x[1], reverse=True)
    # And save
    # Note that here we take into account only
    # novel molecules, i.e. those that are not
    # in ChEMBL (by using their web tool)
    n_top = int(config['BEAM']['n_top'])
    final_mols = []
    final_smiles = []
    for x in sorted_score:
        smiles = x[0]
        if is_novel_chembl(smiles):
            final_smiles.append(smiles)
            final_mols.append(Chem.MolFromSmiles(smiles))
            if len(final_mols)==n_top: 
                break
    
    if final_mols:
        top = Draw.MolsToGridImage(final_mols,
                                   molsPerRow=5,
                                   subImgSize=(400,220),
                                   legends=[str(i+1) for i,x in enumerate(final_mols)])
        top.save(f'{save_path}final_ranking.png')
    else:
        warnings.warn('No valid molecule, so could not output an image')
    
    hp.write_in_file(f'{save_path}final_smiles.txt', final_smiles)
    
    if verbose: print('FINAL RANKING DONE')
    ####################################
    