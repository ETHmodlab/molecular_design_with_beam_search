# Copyright (c) 2021 ETH Zurich

import sys
import argparse
import configparser
import time
import numpy as np

from rdkit.Chem import AllChem
from rdkit import Chem
from sklearn import manifold
import matplotlib.pylab as plt

sys.path.append('../src/')
from python import helper as hp

parser = argparse.ArgumentParser(description='Do t-SNE plot')
parser.add_argument('-c','--configfile', type=str, help='Path to config file', required=True)
parser.add_argument('-v','--verbose', type=bool, help='Verbose', required=True)

def getMorgan(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
    return fp

def get_mol(smi):
    return Chem.MolFromSmiles(smi)

def get_proj(data, n_components=2, perplexity=5):
    """ Function to compute the tSNE"""  
    
    tsne = manifold.TSNE(n_components=n_components, 
                         init='random',
                         random_state=0, 
                         perplexity=perplexity)
    Y = tsne.fit_transform(data)
    
    return Y

def get_canon(smi):
    """
    Input: one SMILES.
    Return: the canonical form
    of the SMILES.
    """
    mol = get_mol(smi)
    if mol is not None: 
        can = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return can
    else:
        return None

def combined_plot(proj, group_idx,
                  m_data="o", s_data=80, alpha_data=0.85,
                  linewidth='0.00', legend=False, shownumber=False):

    fig, ax = plt.subplots(figsize=(10, 10))
    
    margin = 40
    plt.xlim([np.min(proj[:,0])-margin, np.max(proj[:,0])+margin])
    plt.ylim([np.min(proj[:,1])-margin, np.max(proj[:,1])+margin])
    
    labelsize = 16
    plt.xlabel('t-SNE coordinate 1', fontsize=labelsize)
    plt.ylabel('t-SNE coordinate 2', fontsize=labelsize)
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    for name_data,all_idx in group_idx.items():
        plt.scatter(proj[all_idx, 0], proj[all_idx, 1], 
                    lw=0, c=COLOR_PAL_CB[name_data], label=name_data, alpha=alpha_data, s=s_data,
                    marker=m_data, linewidth=linewidth, edgecolors='k')
    
    if legend:
        leg = plt.legend(prop={'size': labelsize}, loc='upper right', markerscale=1.00, frameon=True)
        leg.get_frame().set_alpha(0.9)
        
    return fig
        
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
    
    if verbose: print('\nSTART t-SNE PLOT')
    ####################################
    
    
    
    
    ####################################
    # Path to the data and to save the plots
    dir_exp = f'results/{exp_name}/final_ranking/'
    ####################################
    
    
    
    
    ####################################
    # Get the background data
    path_bground_data = str(config['PLOT']['path_bground_data'])
    
    bground_data = hp.read_with_pd(path_bground_data)
    bground_data = [get_canon(smi) for smi in bground_data]
    
    # the fine-tuning set used
    dir_data = str(config['DATA']['dir_data'])
    name_data = str(config['DATA']['name_data'])
    ft_data = hp.read_with_pd(f'{dir_data}{name_data}')
    
    # top ranked beam search molecules
    final_smi_path = f'{dir_exp}final_smiles.txt'
    try:
        beam_smi = hp.read_with_pd(final_smi_path)
    except:
        print(f'No t-SNE plot - could not find, or empty file at: {final_smi_path}')
        sys.exit(1)
    ####################################
    
    

        
    ####################################
    # do the tSNE
    all_smi = bground_data + ft_data + beam_smi
    all_fp = np.empty([len(all_smi), 1024])
    
    group_idx = {'Background set': [],
                 'Fine-tuning set': [],
                 'Beam designs': []}
    
    for i,smi in enumerate(all_smi):
        mol = get_mol(smi)
        fp = getMorgan(mol)
        arr = np.array(fp)
        all_fp[i] = arr
        
        if smi in beam_smi:
            group_idx['Beam designs'].append(i)
        elif smi in ft_data:
            group_idx['Fine-tuning set'].append(i)
        elif smi in bground_data:
            group_idx['Background set'].append(i)
            
    proj = get_proj(all_fp)
    ####################################
    
    
    
    ####################################
    # do the plot
    COLOR_PAL_CB = {'Background set': 'lightgrey',
                    'Fine-tuning set': '#6666FF',
                    'Beam designs': '#ffa500'}
    
    fig = combined_plot(proj, 
                        group_idx,
                        legend=True)
    fig.savefig(f'{dir_exp}tSNE.png')
        
    end = time.time()
    if verbose: print(f't-SNE PLOT DONE in {end - start:.05} seconds')
    ####################################
        