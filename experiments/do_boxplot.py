# Copyright (c) 2021 ETH Zurich

import sys
import argparse
import configparser
import time
import pandas as pd
import numpy as np
import seaborn as sns
import warnings

from rdkit import DataStructs
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
sys.path.append('../src/python/scaffold_hopping_whales/code/')
import ChemTools as tools
import do_whales

sys.path.append('../src/')
from python import helper as hp
from python import helper_chem as hp_chem

parser = argparse.ArgumentParser(description='Do boxplot')
parser.add_argument('-c','--configfile', type=str, help='Path to config file', required=True)
parser.add_argument('-v','--verbose', type=bool, help='Verbose', required=True)


def get_whales(l_smi):
    """
    Extract WHALES descriptors
    on a list of SMILES
    """
    results = np.zeros([len(l_smi), 33])
    except_idx = []
    for i,smi in enumerate(l_smi):
        result = np.zeros([1, 33])
        template = Chem.MolFromSmiles(smi)
        try:   
            template_prepared, _ = tools.prepare_mol(template)
            whales_template, lab = do_whales.whales_from_mol(template_prepared)
            _whales = np.expand_dims(whales_template, axis=0)
            results[i] = _whales
        except:
            except_idx.append(i)
    
    return results, except_idx

def get_Tanimoto_distance(morgan_a, morgan_b):
    """
    Get the Tanimoto distance between
    two Morgan fingerprints.
    Note: we compute the distance, 
    i.e. 1-similarity
    """
    return 1-DataStructs.FingerprintSimilarity(morgan_a, morgan_b, metric=DataStructs.TanimotoSimilarity)

def get_closest_distance(l_morgan, l_morgan_2):
    results = []
    for a in l_morgan:
        current_best = 1
        for b in l_morgan_2:
            if b!=a:
                distance = get_Tanimoto_distance(a, b)
                if distance<current_best:
                    current_best = distance
        results.append(current_best)
       
    return results

def get_euclidean_distance(l_scaled_w, l_scaled_w_2):
    results = []
    max_dist = -1
    for a in l_scaled_w:
        current_best = None
        for b in l_scaled_w_2:
            if not (a==b).all():
                distance = np.linalg.norm(a-b)
                if not current_best:
                    current_best = distance
                else:
                    if distance<current_best:
                        current_best = distance
                if distance>max_dist:
                    max_dist = distance
        results.append(current_best)
       
    return results, max_dist

def do_boxplot(df, savepath):
    sns.set(rc={'figure.figsize':(8,8)})
    sns.set_style("white")
    sns.color_palette("pastel")
    
    b = sns.boxplot(x="molecule set", hue="Descriptors", y="Distance to background set", 
                    data=df,
                    width=0.35, showfliers=False, linewidth=2.0)
    
    fontsize = 16
    labelsize = 14
    b.set_ylim(0,1)
    b.set_xlabel("molecule set", fontsize=fontsize)
    b.set_ylabel("Distance to background set", fontsize=fontsize)
    b.tick_params(labelsize=labelsize)
    b.legend(fontsize=fontsize) 
    
    sns.despine(top=True, right=True)
    b.figure.savefig(f'{savepath}boxplot.png')
    
    
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
    
    n_core = int(config['PLOT']['n_core'])
    n_top = int(config['BEAM']['n_top'])
    
    if verbose: print('\nSTART BOXPLOT')
    ####################################
    
    
    
    
    ####################################
    # Path to the data and to save the plots
    dir_exp = f'results/{exp_name}/final_ranking/'
    ####################################
    
    
    
    
    ####################################
    # Get the background data
    path_bground_data = str(config['PLOT']['path_bground_data'])
    bground_data = hp.read_with_pd(path_bground_data)
    
    # top ranked beam search molecules
    final_smi_path = f'{dir_exp}final_smiles.txt'
    try:
        beam_smi = hp.read_with_pd(final_smi_path)
    except:
        print(f'No boxplot - could not find, or empty file at: {final_smi_path}')
        sys.exit(1)
    ####################################
    
        
               
        
    ####################################
    # Morgan similarity
    if verbose:
        print('Computing Morgan')
    bground_mols = [hp_chem.get_mol(smi) for smi in bground_data]
    bground_morgan = [hp_chem.get_Morgan(mol) for mol in bground_mols]
    bground_morgan_closest_dist = get_closest_distance(bground_morgan, bground_morgan)
    
    beam_mols = [hp_chem.get_mol(smi) for smi in beam_smi]
    beam_morgan = [hp_chem.get_Morgan(mol) for mol in beam_mols]
    beam_morgan_closest_dist = get_closest_distance(beam_morgan, bground_morgan)
    
    core_morgan_closest_dist = get_closest_distance(beam_morgan[:n_core], bground_morgan)
    
    d_data_m = {'Background set': bground_morgan_closest_dist,
                f'Beam (top {n_top})': beam_morgan_closest_dist,
                f'Beam (top {n_core})': core_morgan_closest_dist}
    ####################################
    
    
    
    
    ####################################
    # WAHLES similarity
    if verbose: print('Computing WHALES')
    bground_w, except_bground_w = get_whales(bground_data)
    if except_bground_w:
        warnings.warn(f'WHALES could not be computed for {len(except_bground_w)} bground data')
    beam_w, except_beam_w = get_whales(beam_smi)
    if except_beam_w:
        warnings.warn(f'WHALES could not be computed for {len(except_beam_w)} beam data')
        
    # we merged them to compute the mean and std,
    # as WHALES need to be scaled
    merged = np.concatenate([bground_w, beam_w])
    stds = np.std(merged, axis=0)
    means = np.mean(merged, axis=0)
    
    scaled_bground_w = (bground_w - means) / stds
    scaled_beam_w = (beam_w - means) / stds
    
    # we get the euclidean distance
    bground_whales_closest_dist, max_d_bg = get_euclidean_distance(scaled_bground_w, scaled_bground_w)
    beams_whales_closest_dist, max_d_b = get_euclidean_distance(scaled_beam_w, scaled_bground_w)
    core_whales_closest_dist, max_d_c = get_euclidean_distance(scaled_beam_w[:n_core], scaled_bground_w)
    # we get the max distance among all
    # to scale values between 0 and 1
    max_dist = np.max([max_d_bg, max_d_b, max_d_c])
    
    # and scale it between 0 and 1 to compare
    # on the same scale on the boxplot
    bground_whales_closest_dist = [x/max_dist for x in bground_whales_closest_dist]
    beams_whales_closest_dist = [x/max_dist for x in beams_whales_closest_dist]
    core_whales_closest_dist = [x/max_dist for x in core_whales_closest_dist]
    
    d_data_w = {'Background set': bground_whales_closest_dist,
                f'Beam (top {n_top})': beams_whales_closest_dist,
                f'Beam (top {n_core})': core_whales_closest_dist}
    ####################################
    
    
    
    
    ####################################
    # let's add everything together in a 
    # df to do a merged boxplot
    if verbose:
        print('Merging data for the plot')
    col_ori = []
    col_des = []
    col_val = []
    
    for k,v in d_data_m.items():
        n_data = len(v)
        col_ori.extend([k]*n_data)
        col_des.extend(['Morgan']*n_data)
        col_val.extend(v)
    
    for k,v in d_data_w.items():
        n_data = len(v)
        col_ori.extend([k]*n_data)
        col_des.extend(['WHALES']*n_data)
        col_val.extend(v)
        
    d_data = {'molecule set': col_ori, 'Descriptors': col_des, 'Distance to background set': col_val}  
    df = pd.DataFrame(d_data)
    ####################################
    
    
    
    
    ####################################
    # do the plot
    do_boxplot(df, dir_exp)
    
    end = time.time()
    if verbose: print(f'BOXPLOT DONE in {end - start:.05} seconds')
    ####################################
        