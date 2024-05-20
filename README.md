# Martini_IDP
Improved Martini force field for IDP

## Introduction:   
1. Martini force field
2. Improved Martini force field
3. However,

## Important tools:   
1. [Martini_OpenMM](https://github.com/maccallumlab/martini_openmm): run Martini simulation in OpenMM   
2. [Martinize2](https://github.com/marrink-lab/vermouth-martinize): automatically convert all-atom to Martini model   
3. [Insane](https://github.com/Tsjerk/Insane): build Martini lipid bilayers ()   

## Tutorials:   
1. [Martini3 protein solution](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2): build protein model in Martini 3   
2. [Martini3 membranes](https://www.sciencedirect.com/science/article/pii/S0076687924000946?via%3Dihub#bib14): build complex membranes with Martini 3 
   
## Build protein-bilayer binding system   
   
example: KR8 + PC_PG bilayer   
### 1. build Martini protein model using Martinize2   
a. installation of Martinize2: ```pip install vermouth```   
b. ```martinize2 -f kr8_at.pdb -o kr8.itp -x kr8_cg.pdb -ff martini3001```   
   
### 2. build bilayer using Insane   
a. installation of Insane: ```pip install insane```   
b. build a 25*25 nm2 bilayer containing 70% POPC and 30% POPG, no salt  
   ```insane -o bilayer.gro -p topol.top -x 25 -y 25 -z 25 -l POPC:7 -l POPG:3 -sol W```   
c. add the following itp files into the topol.top   
    ```#include "martini_v3.0.0.itp"   
        #include "martini_v3.0.0_phospholipids_v1.itp   
        #include "martini_v3.0.0_solvents_v1.itp"   
        #include "martini_v3.0.0_ions_v1.itp"   
        #include "kr8.itp"   
    ```   
c. add salt using gmx genion:   
   ```gmx grompp -f em.mdp -c bilayer.gro -o ions.tpr -p topol.top -maxwarn 1```   
   ```gmx genion -s ions.tpr -p topol.top -o wions.gro -pname NA -nname CL -neutral -conc 0.035```   

### 3. insert protein into the bilayer system   
a. insert 10 kr8 through replace water (W)   
```gmx insert-molecules -f wions.gro -ci kr8_md.gro -nmol 10 -o system.gro -replace W```   
b. update topol.top file   
   get the new number of water: ```grep -c 'W' system.gro```, and change the number of W in *.top file   
   in gromacs 2024, gmx insert-molecules will print out the number of replaced water   
   add the molecule of proteins at the last line: 'KR8      10'   