# Martini_IDP
Improved Martini force field for IDP

## Introduction:   
1. Martini force field
2. Improved Martini force field
3. However,

## Important tools:   
1. [Martini_OpenMM](https://github.com/maccallumlab/martini_openmm): run Martini simulation in OpenMM   
2. [Martinize2](https://github.com/marrink-lab/vermouth-martinize): automatically convert all-atom to Martini model   
3. [Insane](https://github.com/Tsjerk/Insane): build Martini lipid bilayers   

## Tutorials:   
1. [Martini3 protein solution](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2): build protein model in Martini 3   
2. [Martini3 membranes](https://www.sciencedirect.com/science/article/pii/S0076687924000946?via%3Dihub#bib14): build complex membranes with Martini 3 
   
## Build protein-bilayer binding system   
example: KR8 + PC_PG bilayer   
### 1. build Martini protein model using Martinize2   
a. installation of Martinize2: ```pip install vermouth```   
b. ```martinize2 -f kr8_at.pdb -o kr8.itp -x kr8_cg.pdb -ff martini3001```