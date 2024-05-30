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
&emsp;pip install is problometic, use ```python setup.py install```   
4. [PeptideBuilder](https://github.com/clauswilke/PeptideBuilder)   
&emsp;```pip install PeptideBuilder```   

## Tutorials:   
1. [Martini3 protein solution](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2): build protein model in Martini 3   
2. [Martini3 membranes](https://www.sciencedirect.com/science/article/pii/S0076687924000946?via%3Dihub#bib14): build complex membranes with Martini 3 

## Scripts:
1. **scale.py:** scaling the whole martini itp file   
**usage:**```python scale.py scaling_factor```
2. **rename_itp.py:** rename the atomtype name in protein itp file   
**usage:**```python rename_itp.py itp_file```   
3. **run_martini.py:** run Martini simulation in OpenMM   
4. **cg_bond.tcl:** display bonds for Martini model in VMD   
&emsp;```source cg_bond.tcl```   
&emsp;```cg_bonds -topoltype "martini" -top topol.top```   
5. **build_peptide.py:** build all-atom pdb file of IDP   
**usage:** ```python build_peptide.py seq_file```   
**seq_file:** ```name sequence```   

## Examples:
## I. Simulation of protein solution
example: KR8 in 0.035 NaCl solution   
### 1. build Martini protein model
for disordered peptide:   
&emsp;```martinize2 -f kr8_at.pdb -o topol.top -x kr8_cg.pdb -ff martini3001```   
for folded protein:   
&emsp;```martinize2 -f kr8_at.pdb -o topol.top -x kr8_cg.pdb -dssp -ff martini3001 -elastic -ef 700.0 -el 0.5 -eu 0.9 -ea 0 -ep 0 -scfix -cys auto```   
&emsp;dssp library in mdtraj will be used to calculate secondary structure information   
**Note:** martinize2 will also create the molecule_0.itp as default   
**Note:** you should also revise the itp file included in topol.top file   

### 2. build protein solution
a. create box   
&emsp;```gmx editconf -f kr8_cg.pdb -o box.gro -bt cubic -d 1.0```   

b. add water   
&emsp;```gmx solvate -cp box.gro -cs water.gro -radius 0.21 -o boxw.gro -p topol.top```   

c. neutralize and add salt   
&emsp;```gmx grompp -f em.mdp -c boxw.gro -o ions.tpr -p topol.top -maxwarn 1```   
&emsp;```gmx genion -s ions.tpr -p topol.top -o system.gro -pname NA -nname CL -neutral -conc 0.035```   

d. run simulation   
&emsp;**in Gromacs:**
>d1. minimization   
&emsp;```gmx grompp -f em.mdp -c system.gro -p topol.top -o em```   
&emsp;```gmx mdrun -deffem em```   
d2. production   
&emsp;```gmx grompp -f md.mdp -c em.gro -p topol.top -o md```   
&emsp;```gmx mdrun -deffnm md```   

&emsp;**in OpenMM:**
>&emsp;```python run_martini.py```   
&emsp;details can be found below

## II. Build protein-bilayer binding system   
example: KR8 + PC_PG bilayer   
### 1. build Martini protein model using Martinize2   
a. installation of Martinize2: ```pip install vermouth```   
b. convert all-atom model to martini model:   
&emsp;```martinize2 -f kr8_at.pdb -o kr8.itp -x kr8_cg.pdb -ff martini3001```   
   
### 2. build bilayer using Insane   
a. installation of Insane: ```python setup.py install```   
b. build a 25*25 nm2 bilayer containing 70% POPC and 30% POPG, no salt  
&emsp;```insane -o bilayer.gro -p topol.top -x 25 -y 25 -z 25 -l POPC:7 -l POPG:3 -sol W -ff M3```   
c. add the following itp files into the topol.top   
&emsp;```#include "martini_v3.0.0.itp"```   
&emsp;```#include "martini_v3.0.0_phospholipids_v1.itp```   
&emsp;```#include "martini_v3.0.0_solvents_v1.itp"```   
&emsp;```#include "martini_v3.0.0_ions_v1.itp"```   
&emsp;```#include "kr8.itp"```   

### 3. insert protein into the bilayer system   
a. insert 10 kr8 through replace water (W)   
&emsp;```gmx insert-molecules -f bilayer.gro -ci kr8_md.gro -nmol 10 -o bilayer_protein.gro -replace W```   
b. update topol.top file   
&emsp;Get the new number of water: ```grep -c 'W' bilayer_protein.gro```, and change the number of W in *.top file.   
&emsp;In gromacs 2024, "gmx insert-molecules" will print out the number of replaced water   
&emsp;Add the molecule of proteins at the last line: 'KR8      10'   

### 4. add salt
&emsp;add salt using gmx genion:   
&emsp;```gmx grompp -f em.mdp -c bilayer_protein.gro -o ions.tpr -p topol.top -maxwarn 1```   
&emsp;```gmx genion -s ions.tpr -p topol.top -o system.gro -pname NA -nname CL -neutral -conc 0.035```   


## Run Martini in OpenMM   
[martini_openmm](https://github.com/maccallumlab/martini_openmm) is designed for CG simulation in the OpenMM packages.    
To run martini simulation in OpenMM, use [run_martini.py](script/run_martini.py)   
### 1. define simulation system   
&emsp;modify the file name and simulation setups at the begging   
&emsp;for different system, friction constant of LangevinIntegrator should be modified to get the same result with gromacs   
### 2. barostat   
&emsp;modify the barostat in "NPT" section for bilayer or non-bilayer system   
### 3. run simulation
&emsp;```python run_martini.py```   
