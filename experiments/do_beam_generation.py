# Copyright (c) 2021 ETH Zurich

import os, sys
import argparse
import configparser
import time

import numpy as np
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
from rdkit.Chem import Draw

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from keras.models import load_model

sys.path.append('../src/')
from python import helper as hp
from python import fixed_parameters as FP

parser = argparse.ArgumentParser(description='Run beam search')
parser.add_argument('-c','--configfile', type=str, help='Path to config file', required=True)
parser.add_argument('-v','--verbose', type=bool, help='Verbose', required=True)


def int_to_smile(array, indices_token, pad_char):
    """ 
    From an array of int, return a list of 
    molecules in string smile format
    Note: remove the padding char
    """
    all_mols = []
    for seq in array:
        new_mol = [indices_token[str(int(x))] for x in seq]
        all_mols.append(''.join(new_mol).replace(pad_char, ''))
    return all_mols

def one_hot_encode(token_lists, n_chars):
    output = np.zeros((1, len(token_lists), n_chars))
    for i, token in enumerate(token_lists):
        output[0, i, int(token)] = 1
    return output

def save_smiles(candidates, scores, indices_token, start_char, pad_char, end_char, save_path, name_file):
    """
    Save the valid SMILES, along with 
    their score and a picture representation.
    """
    
    all_smi = []
    all_mols = [] # rdkit format
    all_scores = []
    
    # Check SMILES validity
    for x,s in zip(candidates, scores):
        # we have to do one more loop because x is a np array 
        # of dimensions bz, len_smiles, vocab. as the bz is one,
        # which is needed for the keras model, we have to do
        # one more loop to extract the SMILES
        for y in x:
            ints = [indices_token[str(np.argmax(w))] for w in y]
            smiles = ''.join(ints).replace(start_char,'').replace(pad_char,'').replace(end_char,'')
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None: 
                cans = Chem.MolToSmiles(mol, canonical=True)
                if len(cans)>1:
                    all_smi.append(cans)
                    all_mols.append(mol)
                    all_scores.append(s)
    
    hp.write_in_file(f'{save_path}{name_file}_SMILES.txt', all_smi)
    d_smi_to_score = dict(zip(all_smi, all_scores))
    hp.save_obj(d_smi_to_score, f'{save_path}{name_file}_smi_to_score.pkl')
    
    # Draw molecules
    if all_smi:
        top = Draw.MolsToGridImage(all_mols,
                                   molsPerRow=5,
                                   subImgSize=(400,220),
                                   legends=[str(i+1) for i,x in enumerate(all_mols)])
        
        os.makedirs(f'{save_path}img/', exist_ok=True)
        top.save(f'{save_path}img/{name_file}_top.png')

def beam_search_decoder(k, model, vocab_size, max_len, 
                        indices_token, token_indices, name_file, 
                        start_string, pad_char, end_char, save_path, verbose):
    
    seed_token = [token_indices[x] for x in start_string]
    max_len = max_len+1 # to account for the start char
    
    # candidates is a matrix of len(seed_token) 
    # one hot encoded with k elements
    X = one_hot_encode(seed_token, vocab_size)
    candidates = [X]
    scores = [1]*k
    
    for j in range(max_len):
        current_candidates = []
        current_scores = []
        
        for i,x in enumerate(candidates):
            preds = model.predict(x, verbose=0)[0]
            preds = preds[-1, :]
            preds = np.asarray(preds).astype('float64')
            preds = np.log(preds)
            
            # Argsort in descending order, only for the k top
            idx_preds_sorted = np.argsort(preds)[::-1][:k]
            preds_sorted = preds[idx_preds_sorted]
            
            for idx_pred in idx_preds_sorted:
                vec = one_hot_encode([idx_pred], vocab_size)
                new_seq = np.concatenate((x, vec), axis=1)
                current_candidates.append(new_seq)
            
            current_scores.extend([a+b for a,b in zip(scores,preds_sorted.tolist())])
            
        # Find the k best candidates from the scores
        idx_current_best = np.argsort(current_scores)[::-1][:k]
        candidates = [x for i,x in enumerate(current_candidates) if i in idx_current_best]
        scores = [x for i,x in enumerate(current_scores) if i in idx_current_best]           
                
    # save results
    save_smiles(candidates, scores, indices_token, start_char, pad_char, end_char, save_path, name_file)
 
        
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
    
    # get back the experiment parameters
    max_len = int(config['PROCESSING']['max_len'])
    width = int(config['BEAM']['width'])
    
    if verbose: print('\nSTART SAMPLING')
    ####################################
    
    
    
    
    ####################################
    # Path to save the sampled SMILES
    save_path = f'results/{exp_name}/sampling/'
    os.makedirs(save_path, exist_ok=True)
    ####################################
    

    
    
    ####################################
    # Generator parameters
    pad_char = FP.PROCESSING_FIXED['pad_char']
    start_char = FP.PROCESSING_FIXED['start_char']
    end_char = FP.PROCESSING_FIXED['end_char']
    indices_token = FP.INDICES_TOKEN
    token_indices = FP.TOKEN_INDICES
    vocab_size = len(indices_token)
    ####################################
    
    
    
    
    ####################################
    # start iterating over the files
    path_models = f'results/{exp_name}/models/'
    for filename in os.listdir(path_models):
        if filename.endswith('.h5'):
            print(f'--> Currently sampling for model: {filename}')
            path_model = os.path.join(path_models, filename)
            model = load_model(path_model)
            # note: name_file represents the epoch number
            name_file = filename.replace('.h5', '').split('_')[0]
            
            beam_search_decoder(width, model, vocab_size, max_len,
                                indices_token, token_indices, name_file, 
                                start_char, pad_char, end_char, save_path, verbose)
        
            
    end = time.time()
    if verbose: print(f'SAMPLING DONE in {end - start:.05} seconds')
    ####################################
        