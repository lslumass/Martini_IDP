import sys

'''
This script is used for renaming the atom type name in itp file.
Usage: python rename_itp.py itp_file
Date: May 20, 2024
Authour: Shanlong Li@UMass
'''
itp_file = sys.argv[1]

atom_start = None
atom_end = None

with open(itp_file, 'r') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if atom_start is None and line.startswith('[ atoms ]'):
        atom_start = i
        
        for j in range(i+1, len(lines)):
            if lines[j].strip() == '':
                atom_end = j-1
                break
        break

print(atom_start, atom_end)
with open(itp_file[:-4]+'_.itp', 'w') as itp:
    for i, line in enumerate(lines):
        if i <= atom_start or i > atom_end:
            print(line[:-1], file=itp)
        else:
            l = line.split()
            print(l[0].rjust(3)+' '+(l[1]+'_').ljust(5)+' '+line[9:-1], file=itp)

print('Finished!')