import mdtools as mdt 
from mdmodel import mdmodel as mm 
import os,glob

mol = mm()
files = glob.glob('*.xyz')
# file = files[0]
for file in files:
    mol.read_xyz(file)
    # add mass charge
    mol.mass = [[1,12],[2,1]]
    mol.data['charge'] = 0
    mol.natomtype = 2
    mol.data = mdt.transition(mol.data, recenter=True)
    fn = file.split('.')[0]
    mol.write_lmp_data("{}.data".format(fn),'charge',velocity=False)