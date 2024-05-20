import sys

'''
This script is used for scaling martini force field.
Usage: python scale.py scaling_factor
Date: May 20, 2024
Authour: Shanlong Li@UMass
'''
scaling_factor = float(sys.argv[1])

with open('martini_v3.0.0.itp', 'r') as f:
    lines = f.readlines()

atoms = lines[13:857]
nonbs = lines[859:]

with open('martini_v3.0.0_scale.itp', 'w') as f:
    print(';;;; Rebalance Martini 3 Force field for IDPs through scaling protein-protein interaction', file=f)
    print(';;;; scaled atomtypes were labeled with _', file=f)
    print(';;;; created by scale.py\n', file=f)
    print('[ defaults ]\n1 2 ; sigma-epsilon format of LJ parameters\n\n[ atomtypes ]', file=f)

    for atom in atoms:
        print(atom[:-1], file=f)
        a = atom.split()
        idx = len(a[0])
        print(atom[:idx]+'_ '+atom[idx:-1], file=f)
    
    print('\n[ nonbond_params ]', file=f)
    for nonb in nonbs:
        l = nonb.split()
        print(l[0].rjust(6)+'  '+l[1].rjust(6)+'  1 '+l[3]+'    '+l[4], file=f)
        if l[0] == l[1]:
            print(l[0].rjust(6)+'  '+(l[1]+'_').rjust(6)+'  1 '+l[3]+'    '+l[4], file=f)
        else:
            print(l[0].rjust(6)+'  '+(l[1]+'_').rjust(6)+'  1 '+l[3]+'    '+l[4], file=f)
            print((l[0]+'_').rjust(6)+'  '+l[1].rjust(6)+'  1 '+l[3]+'    '+l[4], file=f)
        print((l[0]+'_').rjust(6)+'  '+(l[1]+'_').rjust(6)+'  1 '+l[3]+'    {:.6e}'.format(float(l[4])*scaling_factor), file=f)

print("Finished!")