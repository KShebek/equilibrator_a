# equilibrator_a

## ‼️ IMPORTANT ‼️
### This repo has been incoporated fully into eQuilibrator in the following [repo](https://gitlab.com/KShebek/equilibrator-assets), which is a fork of equilibrator_assets. This code will eventually be merged into the equilibrator master branch. 

### While this repo works, it is advised to use the forked branch linked above.

## Incorporation into eQuilibrator
This repository started out as a way to generate the thermodynamic compounds of biochemical species given a SMILES string. Elad Noor has since incorporated the work I did here into [equilibrator_assets](https://gitlab.com/equilibrator/equilibrator-assets). With [examples of how to use here](https://gitlab.com/equilibrator/equilibrator-assets/-/blob/master/notebooks/generate_compound.ipynb). [Further work](https://gitlab.com/KShebek/equilibrator-assets) done here to incorporate a local database have been integrated into the eQuilibrator environment as well.

## Overview
This repository builds upon arbitrary SMILES support for [equilibrator](https://gitlab.com/equilibrator), the biochemical thermodynamic calculator and incorporates the results into a database for later retrieval. 
This work was the initial place for calculating thermodynamic properties of arbitrary SMILES, but since has been incorporated into eQuilibrator. This repo now highlights how to maintain a local database of compounds for use with the already generated component cache made by Noor hosted on quilt.

## Local component_cache
This code and Noor's code generated thermodynamic properties as follows:
  1. A SMILES string is supplied.
  2. Compound information (pKa, microspecies, etc) is obtained from an existing database or calculated
  3. The compound information is used to calculte various ∆G values. 
  
 What this repo provides is a method to save the results from these calculations into a .sqlite database for later recall.

