# Beam search for automated design and scoring of novel ROR ligands with machine intelligence


## Table of Contents
1. [Description](#Description)
2. [Requirements](#Requirements)
3. [How to run an example experiment](#Run_examples)
5. [How to run the experiments of the paper](#Run_paper)
6. [How to run an experiment on your own data](#Run_own)
7. [Code organisation](#organisation)
8. [Note](#Note)
9. [Acknowledgements](#Acknowledgements)
10. [How to cite this work](#Cite)
11. [License](#license)
12. [Address](#Address)


### Description<a name="Description"></a>

This is the supporting code for the paper «Beam search for automated design and scoring of novel ROR ligands with machine intelligence». This code allows you to replicate the experiments of the paper as well as running our method on your own set of molecules.

Access on the journal webpage (TODO: add the link when officially published)   

[Preprint version (not up to date with the published version)](https://chemrxiv.org/articles/preprint/Beam_Search_Sampling_for_Molecular_Design_and_Intrinsic_Prioritization_with_Machine_Intelligence/14153408?file=26677325)   


**Abstract of the paper**: Chemical language models enable de novo drug design without the requirement for explicit molecular construction rules. While such models have been applied to generate novel compounds with desired bioactivity, the actual prioritization and selection of the most promising computational designs remains challenging. In this work, we leveraged the probabilities learnt by chemical language models with the beam search algorithm as a model-intrinsic technique for automated molecule design and scoring. Prospective application of this method yielded three novel inverse agonists of retinoic acid receptor-related orphan receptors (RORs). Each design was synthesizable in three reaction steps and presented low-micromolar to nanomolar potency towards ROR&gamma;. This model-intrinsic sampling technique eliminates the strict need for external compound scoring functions, thereby further extending the applicability of generative artificial intelligence to data-driven drug discovery.    

### Requirements<a name="Requirements"></a>

First, you need to [clone the repository](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository):

```
git clone git@github.com:ETHmodlab/molecular_design_with_beam_search.git
```
Then, you can run the following command, which will create a conda virtual environement and install all the needed dependencies (if you don't have conda installed, you can get it first by following the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)):   

```
cd molecular_design_with_beam_search/
sh install.sh
```

This command will also install a git submodule in order to use the [WHALES descriptor](https://github.com/grisoniFr/scaffold_hopping_whales). Once the installation is done, you can activate the conda virtual environement:

```
conda activate molb
```
Please note that you will need to activate this virtual environement every time you want to use this project. 

### How to run an example experiment<a name="Run_examples"></a>

Now, you can try quickly the code with a toy experiment by using our [example configuration file](https://github.com/michael1788/molecular_design_with_beam_search/blob/master/experiments/configfiles/0_parameters_example.ini) (which contains explanation of each parameter). By running it with the following command, you will do a fine-tuning experiment with the four natural products modulators of ROR&gamma; used in this paper, for two epochs only (to be fast), and sample molecules with the beam search:

```
cd experiments/
sh run_morty.sh configfiles/0_parameters_example.ini
```

The results of the analysis can be found in *experiments/results/0_parameters_example/*. There, you can find a picture with the top ranked molecules, a reproduction of the paper's figures and a .txt file with the SMILES of the top ranked molecules. Please note that sampled molecules not found in ChEMBL based on a similarity search with their [webresource client](https://github.com/chembl/chembl_webresource_client) are used in the results.

### How to run the experiments of the paper<a name="Run_paper"></a>

If you want to run the same experiments as in the paper, you can run the following for the fine-tuning on the four natural products:   

```
cd experiments/
sh run_morty.sh configfiles/A_experiment_one_fine_tuning_step.ini
```

Or those two experiments in sequence, as the second part of the experiments (as defined by *B2_experiment_two_fine_tuning_steps.ini*) needs the results of the first part (*B1_experiment_two_fine_tuning_steps.ini*), for the experiment with the two-steps fine-tuning:

```
cd experiments/
sh run_morty.sh configfiles/B1_experiment_two_fine_tuning_steps.ini
```

```
sh run_morty.sh configfiles/B2_experiment_two_fine_tuning_steps.ini
```

Note that the experiment with the fine-tuning on the four natural products is fast, even on a CPU. If you don't have 
a GPU, some patience will be needed, even though we provided the pretrained weights of the chemical language model. 
Moreover, make sure you run *B1_experiment_two_fine_tuning_steps.ini* before *B2_experiment_two_fine_tuning_steps.ini*, as *B2_experiment_two_fine_tuning_steps.ini* uses the model trained in *B1_experiment_two_fine_tuning_steps.ini*.

### How to run an experiment on your own data<a name="Run_own"></a>

To do an experiment on your own set of molecules, you will need to create your own configuration file (the *.ini* file). In this file, you can choose your own parameters for the beam serach and the final ranking, as well as give the path your fine-tuning molecules (a *.txt* file with one SMILES string per line).    
Then, you can just run the following command:

```
sh run_morty.sh configfiles/{your_parameter_file_name}.ini
```

You will find the results of your experiment in *experiments/results/{your_parameter_file_name}/*

### Code organisation <a name="organisation"></a>

The main script (*run_morty.sh*) that allows you to run the full experiment with one command can be used separately. If you wish, for example, to only fine-tune a model on your own data, you can run the following:

```
sh run_training.sh configfiles/{your_parameter_file_name}.ini
```
All specific scripts (to fine-tune, do the plots, etc) can be run in the same way.

### Note <a name="Note"></a>

This work (code and paper) is build on top of our [previous research](https://www.nature.com/articles/s42256-020-0160-y.epdf?author_access_token=kx71VwOu26XWGELCg3BP-NRgN0jAjWel9jnR3ZoTv0MojvyIaQWNqzF7aemIUbYlNUc8tqoGgWco3JoR6d8H9plcxmpko09VfAUvw6-sCHyp8bABy7FhZ89AUc_da9ZU3s4YWQy4gK0meFq2XLhHYA%3D%3D). Notably, if you wish to pretrain a chemical langauge model on your own data—rather than using one of the two available pretrained models here—we recommend you to use the open source code of our previous paper  (https://github.com/ETHmodlab/virtual_libraries).

### Acknowledgements <a name="Acknowledgements"></a>

This research was supported by the Swiss National Science Foundation (grant no. 205321_182176 to Gisbert Schneider), the RETHINK initiative at ETH Zurich and the Novartis Forschungsstiftung (FreeNovationgrant “AI in Drug Discovery” to Gisbert Schneider).
 
### How to cite this work<a name="Cite"></a>

TODO: update with the journal version when published
```
@article{moret2021beam,
  title={Beam search sampling for molecular design and intrinsic prioritization with machine intelligence},
  author={Moret, Michael and Helmst{\"a}dter, Moritz and Grisoni, Francesca and Schneider, Gisbert and Merk, Daniel}
}
```

### License<a name="License"></a>
[MIT License](LICENSE)

### Address<a name="Address"></a>
MODLAB   
ETH Zurich   
Inst. of Pharm. Sciences   
HCI H 413   
Vladimir-​Prelog-Weg 4   
CH-​8093 Zurich   
